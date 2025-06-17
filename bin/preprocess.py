#!/usr/bin/env python3
# Roger Volden

import os
import sys
import mappy
import subprocess
import tempfile
import re

def preprocess(args, tmp_dir, fastq_file):
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    tmp_adapter_dict={}
#    print('\tAligning splints to reads with blat', file=sys.stderr)

    process(args, fastq_file, tmp_dir)

    adapter_set = set()
    with open(align_psl) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            line = line.split('\t')
            read_name, adapter, strand = line[9], line[13], line[8]
            gaps, score = float(line[5]), float(line[0])
            if gaps < 50 and score > 50:
                if read_name not in tmp_adapter_dict:
                    tmp_adapter_dict[read_name]=[[None, 1, None]]
                tmp_adapter_dict[read_name].append([adapter, float(line[0]), strand])
                adapter_set.add(adapter)

    adapter_dict = {} # read_id: [adapter, strand]
    no_splint_reads = 0
    for name, alignments in tmp_adapter_dict.items():
        best = sorted(alignments, key=lambda x: x[1], reverse=True)[0]
        if not best[0]:
            no_splint_reads += 1
            continue
        adapter_set.add(best[0])
        adapter_dict[name] = [best[0], best[2]]
    return adapter_dict, adapter_set, no_splint_reads

def process(args, fastq_file, tmp_dir):
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    
    # Create direct minimap2 alignment with --score-N=0 for degenerate bases
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'
    
    # Load splint sequences
    splint_seqs = {}
    for name, seq, _ in mappy.fastx_read(args.splint_file):
        splint_seqs[name] = seq.upper()
    
    # Load read sequences
    read_seqs = {}
    for name, seq, _ in mappy.fastx_read(fastq_file):
        read_seqs[name] = seq.upper()
    
    # Open output file in PSL-like format
    with open(align_psl, 'w') as psl_out:
        # Process each read against each splint using direct minimap2
        for read_name, read_seq in read_seqs.items():
            for splint_name, splint_seq in splint_seqs.items():
                
                # Create temporary files for minimap2
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as read_file:
                    read_file.write(f">{read_name}\n{read_seq}\n")
                    read_path = read_file.name
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as splint_file:
                    splint_file.write(f">{splint_name}\n{splint_seq}\n")
                    splint_path = splint_file.name
                
                try:
                    # Run minimap2 with --score-N=0 to handle degenerate bases properly
                    # Find minimap2 in conda environment or system PATH
                    minimap2_cmd = 'minimap2'
                    # Check common conda locations
                    conda_paths = [
                        '/Users/gordon/miniforge3/envs/maple/bin/minimap2',
                        '/opt/conda/bin/minimap2',
                        '/usr/local/bin/minimap2'
                    ]
                    for path in conda_paths:
                        if os.path.exists(path):
                            minimap2_cmd = path
                            break
                    
                    cmd = [
                        minimap2_cmd,
                        '--score-N=0',  # Score N bases as 0 (no penalty/bonus)
                        '-ax', 'sr',    # Short read preset
                        '-k', '8',      # Small k-mer for sensitivity
                        '-w', '4',      # Small window
                        '--secondary=yes',  # Report secondary alignments
                        read_path,      # Reference (read sequence)
                        splint_path     # Query (splint sequence)
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    
                    # Parse SAM output and convert to PSL-like format
                    for line in result.stdout.strip().split('\n'):
                        if line.startswith('@') or not line.strip():
                            continue
                        
                        fields = line.split('\t')
                        if len(fields) < 11:
                            continue
                        
                        flag = int(fields[1])
                        if flag & 4:  # Unmapped
                            continue
                        
                        # Parse alignment details
                        pos_start = int(fields[3]) - 1  # Convert to 0-based
                        cigar = fields[5]
                        
                        # Parse CIGAR string
                        cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                        
                        # Calculate alignment metrics
                        ref_length = 0
                        matches = 0
                        gaps = 0
                        
                        for length, op in cigar_ops:
                            length = int(length)
                            if op in 'MDN=X':  # Operations that consume reference
                                ref_length += length
                            if op in 'M=':  # Matches
                                matches += length
                            elif op in 'ID':  # Indels
                                gaps += length
                        
                        pos_end = pos_start + ref_length
                        score = matches  # use matches as score
                        strand = '+' if not (flag & 16) else '-'  # Check reverse complement flag
                        
                        # Write PSL-like line compatible with existing parser
                        # PSL format: matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert,
                        #            tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, 
                        #            tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts
                        psl_line = f"{score}\t0\t0\t0\t0\t{gaps}\t0\t0\t{strand}\t{read_name}\t{len(read_seq)}\t{pos_start}\t{pos_end}\t{splint_name}\t{len(splint_seq)}\t0\t{len(splint_seq)}\t1\t{ref_length}\t{pos_start}\t0\n"
                        psl_out.write(psl_line)
                
                except subprocess.CalledProcessError as e:
                    # Silently continue if minimap2 fails for this read/splint pair
                    pass
                except FileNotFoundError:
                    print(f"Warning: minimap2 not found. Please install minimap2 or add it to PATH.", file=sys.stderr)
                    break
                
                finally:
                    # Clean up temporary files
                    try:
                        os.unlink(read_path)
                        os.unlink(splint_path)
                    except:
                        pass

