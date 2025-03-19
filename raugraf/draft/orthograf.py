import os
import re
import json
from collections import defaultdict
import itertools
import sys
from sys import stderr
from collections import defaultdict
from statistics import mean
from itertools import combinations
from math import ceil, pow


import gffutils


def parse_chromosome_info(chrom_name):
    """Parse sample, haplotype, and chromosome from the chromosome name format ${SAMPLE}_${HAP}_${CHROM}"""
    parts = chrom_name.split('_')
    if len(parts) >= 3:
        chrom = parts[-1]
        hap = parts[-2]
        sample = '_'.join(parts[:-2])  # Handle sample names that might contain underscores
        return sample, hap, chrom
    else:
        raise ValueError("Your contigs aren't named according to the PanSN naming spec, or --pansn-chr is wrong")

def load_orthogroup_mapping(tsv_file, mode="diamond"):
    """Load gene ID to orthogroup ID mapping from TSV file"""
    print(f"Loading orthogroup mapping from {tsv_file}", file=stderr)
    gid2oid = {}
    with open(tsv_file) as file:
        if mode == "diamond":
            for line in file:
                F = line.rstrip().split('\t')
                gid2oid[F[1]] = F[0]
        elif mode == "orthofinder":
            for line in file:
                if line.startswith("Orthogroup\t"):
                    continue
                F = line.rstrip("\r\n").split("\t", 1)
                oid = F[0]
                gids = re.split(", |\t", F[1])
                for gid in gids:
                    if gid:
                        gid2oid[gid] = oid
    print(f"Loaded mapping of {len(gid2oid)} genes to {len(set(gid2oid.values()))} orthogroups", file=stderr)
    return gid2oid

def extract_gene_order_from_gff(gff_file, orthogroup_map, feature_type='gene', pansn_chr="_"):
    """Extract ordered list of orthogroups from a GFF3 file based on gene positions"""
    print(f"Processing GFF file: {gff_file}", file=stderr)
    
    # Create temporary database for GFF parsing
    db_file = f"/dev/shm/{os.path.basename(gff_file)}.db"
    if not os.path.exists(db_file):
        db = gffutils.create_db(gff_file, dbfn=db_file, force=True, merge_strategy='merge',
                               sort_attribute_values=True)
    else:
        db = gffutils.FeatureDB(db_file)
    
    # Dictionary to store chromosome -> list of (position, orthogroup, gene_info) tuples
    # gene_info now includes start and end positions
    chromosome_genes = defaultdict(list)
    missing_genes = 0
    total_genes = 0
    
    # Extract genes and their positions
    for feature in db.features_of_type(feature_type):
        total_genes += 1
        
        # Get gene ID (might be in ID or Name attribute)
        gene_id = feature.id
        if gene_id.startswith(feature_type + ':'):
            gene_id = gene_id[len(feature_type + ':'):]
        
        # Get chromosome/contig name
        chromosome = feature.seqid
        
        # Get gene position (use start position for sorting)
        position = feature.start
        
        # Look up orthogroup
        if gene_id in orthogroup_map:
            orthogroup = orthogroup_map[gene_id]
            # Store additional gene information including start and end positions
            gene_info = {
                'id': gene_id,
                'start': feature.start,
                'end': feature.end,
                'strand': feature.strand,
                'chrom': chromosome,
            }
            chromosome_genes[chromosome].append((position, orthogroup, gene_info))
        else:
            missing_genes += 1
    
    # Sort genes by position for each chromosome
    for chromosome in chromosome_genes:
        chromosome_genes[chromosome].sort(key=lambda x: (x[0], x[1]))  # Sort by position
    
    # Extract final paths (orthogroup IDs and gene positions)
    chromosome_paths = {}
    for chromosome, genes in chromosome_genes.items():
        # Extract sample, haplotype, and chromosome information
        sample, hap, chrom = parse_chromosome_info(chromosome)
        
        # Get orthogroup IDs in order and track gene positions
        orthogroup_path = []
        gene_positions = []
        
        for _, og, gene_info in genes:
            orthogroup_path.append(og)
            gene_positions.append(gene_info)
        
        # Store path with standardized key for grouping later
        chromosome_paths[chromosome] = {
            'sample': sample,
            'haplotype': hap,
            'chromosome': chrom,
            'path': orthogroup_path,
            'gene_positions': gene_positions  # Store gene positions for each orthogroup
        }
    
    print(f"Processed {total_genes} genes, {missing_genes} genes not found in orthogroup mapping", file=stderr)
    print(f"Created paths for {len(chromosome_paths)} chromosomes", file=stderr)
    
    return chromosome_paths

def group_paths_by_chromosome(all_paths):
    """Group paths by chromosome across all samples"""
    chromosome_groups = defaultdict(list)
    
    for _, path_info in all_paths.items():
        chrom = path_info['chromosome']
        chromosome_groups[chrom].append({
            'sample': path_info['sample'],
            'haplotype': path_info['haplotype'],
            'path': path_info['path'],
            'gene_positions': path_info['gene_positions']  # Include gene positions
        })
    
    return chromosome_groups

def find_segment_occurrences(path, segment):
    """Find all start indices where the segment occurs in the path"""
    occurrences = []
    segment_len = len(segment)
    path_len = len(path)
    
    for i in range(path_len - segment_len + 1):
        if path[i:i+segment_len] == segment:
            occurrences.append(i)
    
    return occurrences

def build_segment_index(paths, min_length):
    """
    Build an index of all segments of minimum length and their occurrences in paths.
    
    Args:
        paths: List of paths
        min_length: Minimum segment length to index
        
    Returns:
        Dictionary mapping segments to the list of path indices containing them
    """
    segment_index = defaultdict(list)
    
    # For each path, build a suffix array-like structure
    for path_idx, path in enumerate(paths):
        path_length = len(path)
        
        # Use a sliding window approach with a hash table for indexing
        for start in range(path_length - min_length + 1):
            # Store segments of increasing length starting at this position
            for length in range(min_length, path_length - start + 1):
                segment = path[start:start+length]
                
                # Only add this path_idx once per segment
                if path_idx not in segment_index[segment]:
                    segment_index[segment].append(path_idx)
    
    return segment_index

def filter_contained_segments_with_positions(colinear_segments):
    """
    Filter out segments that are contained within larger segments with the same path occurrences.
    
    Args:
        colinear_segments: Dictionary mapping segments to path indices and positions
        
    Returns:
        Dictionary with only maximally shared segments
    """
    maximally_shared_segments = {}
    sorted_segments = sorted(colinear_segments.keys(), key=len, reverse=True)
    
    for i, segment in enumerate(sorted_segments):
        # Check if this segment is contained in any of the larger segments
        is_contained = False
        segment_paths = set(colinear_segments[segment]['path_indices'])
        
        # Only need to check against larger segments (which come earlier in the sorted list)
        for j in range(i):
            larger_segment = sorted_segments[j]
            larger_segment_paths = set(colinear_segments[larger_segment]['path_indices'])
            
            # If segment is a substring of larger_segment and has the same path occurrences
            if _is_subpath(segment, larger_segment) and segment_paths.issubset(larger_segment_paths):
                is_contained = True
                break
        
        if not is_contained:
            maximally_shared_segments[segment] = colinear_segments[segment]
    
    return maximally_shared_segments

def _is_subpath(smaller, larger):
    """Check if smaller is a contiguous subpath of larger"""
    n, m = len(smaller), len(larger)
    if n > m:
        return False
    
    for i in range(m - n + 1):
        if larger[i:i+n] == smaller:
            return True
    return False

# Modified optimized version to include gene positions
def detect_colinear_segments_optimized(paths, gene_positions_list, min_shared_paths=2, min_length=3):
    """
    Optimized version that reduces memory usage and handles repeating orthogroups efficiently.
    
    Args:
        paths: List of paths, where each path is a list of orthogroup IDs
        gene_positions_list: List of gene position information for each path
        min_shared_paths: Minimum number of paths that must share a segment
        min_length: Minimum length of shared segments to report
        
    Returns:
        Dictionary mapping path segments to path indices and position information
    """
    # Convert paths to tuples for hashability
    tuple_paths = [tuple(path) for path in paths]
    
    # Dictionary to map segments to (path_idx, start_position) pairs
    segment_positions = defaultdict(list)
    colinear_segments = {}
    
    # For each path, record the starting positions of each min_length segment
    for path_idx, path in enumerate(tuple_paths):
        path_length = len(path)
        
        # First pass: index all minimum-length segments
        for start in range(path_length - min_length + 1):
            min_segment = path[start:start+min_length]
            segment_positions[min_segment].append((path_idx, start))
    
    # Process only the segments that appear in enough paths
    for min_segment, positions in segment_positions.items():
        # Group positions by path_idx
        positions_by_path = defaultdict(list)
        for path_idx, start_pos in positions:
            positions_by_path[path_idx].append(start_pos)
        
        # Only process segments that appear in enough different paths
        if len(positions_by_path) >= min_shared_paths:
            # Check if this is a repeating pattern (all elements in min_segment are the same)
            is_repeating = len(set(min_segment)) == 1
            
            # For repeating patterns, pre-process the positions to avoid combinatorial explosion
            if is_repeating:
                # For each path, only keep the first occurrence of the repeating pattern
                # This avoids generating all combinations of repeating positions
                for path_idx in positions_by_path:
                    # Sort positions to ensure we get the earliest valid position
                    positions_by_path[path_idx] = [min(positions_by_path[path_idx])]
            
            # Limit the number of positions per path to prevent explosion
            # even for non-repeating patterns
            max_positions_per_path = 3
            for path_idx in positions_by_path:
                if len(positions_by_path[path_idx]) > max_positions_per_path:
                    # Keep only the first few positions for each path
                    positions_by_path[path_idx] = positions_by_path[path_idx][:max_positions_per_path]
            
            # Select the paths that have enough occurrences for analysis
            valid_path_idxs = [path_idx for path_idx in positions_by_path 
                               if len(positions_by_path[path_idx]) > 0]
            
            # Check if we have enough paths with this segment
            if len(valid_path_idxs) >= min_shared_paths:
                # For each possible combination of paths that meet the minimum count
                for path_idxs in combinations(valid_path_idxs, min_shared_paths):
                    # Get combinations of starting positions
                    pos_combinations = _get_position_combinations(positions_by_path, path_idxs)
                    
                    # Process each position combination
                    for pos_combination in pos_combinations:
                        # For repeating patterns, extend the segment intelligently
                        if is_repeating:
                            max_segment = _find_max_repeating_segment(tuple_paths, path_idxs, 
                                                                      pos_combination, min_segment[0])
                        else:
                            max_segment = _find_max_shared_segment(tuple_paths, path_idxs, 
                                                                  pos_combination, min_length)
                        
                        if max_segment and len(max_segment) >= min_length:
                            # Store segment information
                            if max_segment not in colinear_segments:
                                colinear_segments[max_segment] = {
                                    'path_indices': list(path_idxs),
                                    'positions': []
                                }
                            
                            # Store position information for each path
                            for idx, path_idx in enumerate(path_idxs):
                                start_idx = pos_combination[idx]
                                end_idx = start_idx + len(max_segment) - 1
                                
                                # Add position information if not already present
                                pos_info = {
                                    'path_idx': path_idx,
                                    'start_idx': start_idx,
                                    'end_idx': end_idx,
                                    'start_pos': gene_positions_list[path_idx][start_idx]['start'],
                                    'end_pos': gene_positions_list[path_idx][end_idx]['end'],
                                    'chrom': gene_positions_list[path_idx][end_idx]['chrom'],
                                }
                                
                                # Check if this position info is already recorded
                                already_recorded = False
                                for existing_pos in colinear_segments[max_segment]['positions']:
                                    if (existing_pos['path_idx'] == path_idx and 
                                        existing_pos['start_idx'] == start_idx):
                                        already_recorded = True
                                        break
                                
                                if not already_recorded:
                                    colinear_segments[max_segment]['positions'].append(pos_info)
    
    # Filter out segments contained in larger segments
    return filter_contained_segments_with_positions(colinear_segments)

def _get_position_combinations(positions_by_path, path_idxs):
    """Generate combinations of starting positions across multiple paths with safety limits"""
    position_lists = [positions_by_path[idx] for idx in path_idxs]
    
    # Safety check - if total combinations would be too large, reduce the number of positions
    total_combinations = 1
    for pos_list in position_lists:
        total_combinations *= len(pos_list)
    
    # If combinations would exceed threshold, limit each list further
    max_safe_combinations = 100
    if total_combinations > max_safe_combinations:
        reduction_factor = int(ceil(pow(total_combinations / max_safe_combinations, 1/len(position_lists))))
        position_lists = [lst[:max(1, len(lst) // reduction_factor)] for lst in position_lists]
    
    return itertools.product(*position_lists)

def _find_max_repeating_segment(paths, path_idxs, positions, repeating_element):
    """
    Efficiently find the longest repeating segment starting at given positions.
    Specialized for segments with the same element repeated.
    """
    # Get the paths we're comparing
    relevant_paths = [paths[idx] for idx in path_idxs]
    path_lengths = [len(path) for path in relevant_paths]
    
    # Count consecutive occurrences of the repeating element in each path
    repeating_counts = []
    for path_idx, pos in enumerate(positions):
        count = 0
        path = relevant_paths[path_idx]
        
        # Count consecutive occurrences of the same element
        while pos + count < path_lengths[path_idx] and path[pos + count] == repeating_element:
            count += 1
        
        repeating_counts.append(count)
    
    # The maximum length of the repeating segment is the minimum count across all paths
    max_length = min(repeating_counts)
    
    if max_length > 0:
        # Create the segment of repeated elements
        return tuple([repeating_element] * max_length)
    return None


def _find_max_shared_segment(paths, path_idxs, positions, min_length):
    """Find the longest shared segment starting at given positions in the paths"""
    # Get the paths we're comparing
    relevant_paths = [paths[idx] for idx in path_idxs]
    path_lengths = [len(path) for path in relevant_paths]
    
    # Calculate max possible length
    max_possible_length = min(length - pos for length, pos in zip(path_lengths, positions))
    
    # Early exit when all paths have identical content at the starting positions
    if max_possible_length < min_length:
        return None
    
    # Check if we're starting with a repeating pattern (same OG repeated)
    start_elements = [path[pos] for path, pos in zip(relevant_paths, positions)]
    if all(el == start_elements[0] for el in start_elements):
        # Check for repeating elements (e.g., AAAA pattern)
        repeating_counts = []
        for path_idx, pos in enumerate(positions):
            count = 1
            path = relevant_paths[path_idx]
            element = path[pos]
            # Count consecutive occurrences of the same element
            while pos + count < path_lengths[path_idx] and path[pos + count] == element:
                count += 1
            repeating_counts.append(count)
        
        # If we have repeating elements, limit the extension to the minimum repetition
        if min(repeating_counts) > 1:
            # Limit the max_possible_length to the minimum number of repeats
            consecutive_repeats = min(repeating_counts)
            # Force the segment to break after the repeats to avoid infinite loops
            max_possible_length = min(max_possible_length, consecutive_repeats)
    
    # Find maximum length where all segments match
    max_length = 0
    reference_segment = None
    
    for length in range(min_length, max_possible_length + 1):
        current_segments = [path[pos:pos+length] for path, pos in zip(relevant_paths, positions)]
        # Check if all segments match the first one
        if all(seg == current_segments[0] for seg in current_segments):
            max_length = length
            reference_segment = current_segments[0]
        else:
            break
    
    if max_length >= min_length:
        return reference_segment
    return None

def identify_neighbourhoods(colinear_segments, sample_info, gene_positions_list, just_chrom):
    """
    Identify neighbourhoods (non-syntenic regions) between syntenic blocks.
    
    Args:
        colinear_segments: Dictionary of colinear segments
        sample_info: List of sample names
        gene_positions_list: List of gene position information for each path
        
    Returns:
        List of neighbourhood regions
    """
    # First, organize segments by sample and sort by position
    segments_by_sample = defaultdict(list)
    
    for segment, segment_info in colinear_segments.items():
        for pos_info in segment_info['positions']:
            sample_idx = pos_info['path_idx']
            sample_name = sample_info[sample_idx]
            
            segments_by_sample[sample_idx].append({
                'segment': segment,
                'start_pos': pos_info['start_pos'],
                'end_pos': pos_info['end_pos'],
                'start_idx': pos_info['start_idx'],
                'end_idx': pos_info['end_idx'],
                'chrom': pos_info['chrom'],
                'sample': sample_name
            })
    
    # Sort segments within each sample by position
    for sample_idx in segments_by_sample:
        segments_by_sample[sample_idx].sort(key=lambda x: x['start_pos'])
    
    # Identify neighbourhoods between segments
    neighbourhoods = []
    neighbourhood_id = 1
    
    # Get all unique pairs of adjacent segments across all samples
    segment_pairs = set()
    
    for sample_idx, segments in segments_by_sample.items():
        for i in range(len(segments) - 1):
            # Create a unique identifier for this segment pair
            # Use the segment tuples themselves as identifiers
            left_segment = segments[i]['segment']
            right_segment = segments[i+1]['segment']
            
            # Only add if not already processed
            segment_pairs.add((left_segment, right_segment))
    
    # Process each unique segment pair to create neighbourhoods
    for left_segment, right_segment in segment_pairs:
        neighbourhood = {
            'id': f"{just_chrom}_nh{neighbourhood_id:05d}",
            'left_segment': list(left_segment),
            'right_segment': list(right_segment),
            'regions': []
        }
        
        # Find all instances of this segment pair across samples
        for sample_idx, segments in segments_by_sample.items():
            for i in range(len(segments) - 1):
                if segments[i]['segment'] == left_segment and segments[i+1]['segment'] == right_segment:
                    # This is an instance of our neighbourhood
                    # The neighbourhood extends from the end of the left segment to the start of the right segment
                    region = {
                        'sample': segments[i]['sample'],
                        'chrom': segments[i]['chrom'],
                        'start_pos': segments[i]['end_pos'],  # End of left syntenic block
                        'end_pos': segments[i+1]['start_pos'],  # Start of right syntenic block
                        'sample_idx': sample_idx
                    }
                    neighbourhood['regions'].append(region)
        
        # Only add neighbourhoods that appear in at least one sample
        if neighbourhood['regions']:
            neighbourhoods.append(neighbourhood)
            neighbourhood_id += 1
    
    return neighbourhoods

def orthograf(gff_files, orthogroup_map, min_shared_paths=3, min_segment_length=5, bedfile=None, neighbourhood_file=None, pansn_chr="_"):
    """Main function to process GFF files and analyze orthogroup paths"""
    # Process all GFF files provided on the command line
    all_paths = {}
    for gff_path in gff_files:
        if not os.path.exists(gff_path):
            print(f"Warning: GFF file '{gff_path}' does not exist. Skipping.", file=stderr)
            continue
            
        chromosome_paths = extract_gene_order_from_gff(gff_path, orthogroup_map, pansn_chr=pansn_chr)
        all_paths.update(chromosome_paths)
    
    chromosome_groups = group_paths_by_chromosome(all_paths)
    
    results = {}
    print("\n")
    for chrom, paths_info in chromosome_groups.items():
        if len(paths_info) < min_shared_paths:
            print(f"Skip {chrom} as it has only {len(paths_info)} samples", file=stderr)
            continue

        print(f"Analyzing chromosome {chrom} across {len(paths_info)} samples/haplotypes", file=stderr)
        
        just_paths = [p['path'] for p in paths_info]
        gene_positions_list = [p['gene_positions'] for p in paths_info]
        sample_info = [f"{p['sample']}_{p['haplotype']}" for p in paths_info]
        
        colinear_segments = detect_colinear_segments_optimized(
            just_paths, 
            gene_positions_list,
            min_shared_paths=min_shared_paths, 
            min_length=min_segment_length
        )
        
        neighbourhoods = identify_neighbourhoods(colinear_segments, sample_info, gene_positions_list, just_chrom=chrom)
        
        results[chrom] = {
            'sample_count': len(paths_info),
            'samples': sample_info,
            'colinear_segments': [
                {
                    'segment': list(segment),
                    'length': len(segment),
                    'shared_by': len(segment_info['path_indices']),
                    'samples': [sample_info[idx] for idx in segment_info['path_indices']],
                    'positions': [
                        {
                            'sample': sample_info[pos_info['path_idx']],
                            'start_pos': pos_info['start_pos'],
                            'end_pos': pos_info['end_pos'],
                            'chrom': pos_info['chrom'],
                        }
                        for pos_info in segment_info['positions']
                    ]
                }
                for segment, segment_info in colinear_segments.items()
            ],
            'neighbourhoods': neighbourhoods
        }
        
        # Print summary
        print(f"Found {len(colinear_segments)} co-linear segments shared by at least {min_shared_paths} samples", file=stderr)
        print(f"Found {len(neighbourhoods)} neighbourhoods between syntenic blocks", file=stderr)
        print("\n")
    
    if bedfile:
        for chrom, result in results.items():
            for i, seg in enumerate(sorted(result['colinear_segments'], key=lambda x: mean(y["start_pos"] for y in x["positions"]))):
                print(f"\n# Segment of {seg['length']} orthogroups shared by {seg['shared_by']} samples", file=bedfile)
                print("#", " -> ".join(seg['segment']), file=bedfile)
                for pos in seg['positions']:
                    print(pos['chrom'], pos['start_pos'], pos['end_pos'], sep="\t", file=bedfile)
    
    if neighbourhood_file:
        for chrom, result in results.items():
            for nh in result['neighbourhoods']:
                print(f"\n# Neighbourhood {nh['id']} between syntenic blocks", file=neighbourhood_file)
                print("# Left block:", " -> ".join(nh['left_segment']), file=neighbourhood_file)
                print("# Right block:", " -> ".join(nh['right_segment']), file=neighbourhood_file)
                for region in nh['regions']:
                    print(region['chrom'], region['start_pos'], region['end_pos'], nh['id'], sep="\t", file=neighbourhood_file)

    print("\nAnalysis complete", file=stderr)
    return results

def main(argv=None):
    """Orthograf: generate neighbourhoods where synteny is broken.

    DRAFT TOOL, use with care"""

    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze orthogroup paths from GFF files')
    parser.add_argument('gff_files', nargs='+', help='GFF files to process')
    parser.add_argument('--ogs', help='GeneID-OGID mapping file from either diamond or orthofinder')
    parser.add_argument('--og-format', help='Which format is --ogs?', choices=("diamond", "orthofinder"), default="diamond")
    parser.add_argument('--min-shared', type=int, default=17, help='Minimum number of samples sharing a segment')
    parser.add_argument('--min-length', type=int, default=5, help='Minimum length of co-linear segments in genes')
    parser.add_argument('--out-syntbed', help="Bedfile output for syntenic regions", type=argparse.FileType("wt"))
    parser.add_argument('--out-nhbed', help="Bedfile output for neighbourhood regions", type=argparse.FileType("wt"))
    parser.add_argument('--out-json', help="JSON data dump output", type=argparse.FileType("wt"))
    parser.add_argument('--pansn-chr', default="_", help="Separator between sample, hap, and chrom in PanSN names (see PanSN spec, default is _ NOT #)")
    
    args = parser.parse_args(argv)
    og_mapping = load_orthogroup_mapping(args.ogs, args.og_format)
    res = orthograf(args.gff_files, og_mapping, args.min_shared, args.min_length, 
               bedfile=args.out_syntbed, neighbourhood_file=args.out_nhbed,
               pansn_chr=args.pansn_chr)
    if args.out_json is not None:
        json.dump(res, args.out_json, indent=4)

if __name__ == "__main__":
    main()
