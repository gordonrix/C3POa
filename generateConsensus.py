import os
import sys
import numpy as np
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import mappy as mm
from conk import conk
import gc
import gzip
import time
import subprocess

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

C3POaPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
abpoa='/Users/gordon/miniforge3/envs/maple/src/abPOA-v1.5.3/bin/abpoa'
racon='/Users/gordon/miniforge3/envs/maple/bin/racon'
# BLAT removed - using minimap2/mappy instead

from preprocess import preprocess
from call_peaks import call_peaks
from determine_consensus import determine_consensus

def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x) / base))


def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Makes consensus sequences from R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--reads', '-r', type=str, action='store',
                          help='FASTQ file that contains the long R2C2 reads or a folder containing multiple of these FASTQ files.')
    parser.add_argument('--splint_file', '-s', type=str, action='store',
                          help='Path to the splint FASTA file.')
    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='''Sets the length cutoff for your raw sequences. Anything
                                shorter than the cutoff will be excluded. Defaults to 1000.''')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='''Sets the median distance cutoff for consensus sequences.
                                Anything shorter will be excluded. Defaults to 500.''')
    parser.add_argument('--nosplint', '-ns', action='store_true',
                        help='''When set the first 300 bases of each read are used as splint''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--peakFinderSettings', '-p', action='store', default='23,3,35,2')
    parser.add_argument('--batch_size', '-bs', type=int, default=1000,
                        help='Number of reads to batch together for abpoa processing. Defaults to 1000.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def analyze_reads(args, read, splint, read_adapter, racon, tmp_dir,abpoa):

    peakFinderSettings=args.peakFinderSettings.split(',')
    penalty, iters, window, order = int(peakFinderSettings[0]),int(peakFinderSettings[1]),int(peakFinderSettings[2]),int(peakFinderSettings[3])

    name, seq, qual = read[0], read[1], read[2]
    seq_len = len(seq)
    use=False
    read_consensus = ''
    subs=[]
    
    # Time the conk alignment
    conk_start = time.time()
    scores = conk.conk(splint, seq, penalty)
    conk_time = time.time() - conk_start
    
    # Time the peak detection
    peaks_start = time.time()
    peaks = call_peaks(scores, args.mdistcutoff, iters, window, order)
    peaks_time = time.time() - peaks_start
    if  list(peaks):
        peaks = list(peaks + len(splint) // 2)
        for i in range(len(peaks) - 1, -1, -1):
            if peaks[i] >= seq_len:
                del peaks[i]
        if peaks:

            if args.nosplint:
                if len(peaks)>1:
                    subread_lens = np.diff(peaks)
                    relative_std = np.std(subread_lens)/np.mean(subread_lens)
                    if relative_std < 0.1:
                        use=True
                else:
                    use=True
            else:
                use=True

        # check for outliers in subread length
    if use:
        subreads, qual_subreads, dangling_subreads, qual_dangling_subreads = [], [], [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks)
            subread_lens = [rounding(x, 50) for x in subread_lens]
            median_subread_len = np.median(subread_lens)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2:
                    subreads.append(seq[bounds[0]:bounds[1]])
                    qual_subreads.append(qual[bounds[0]:bounds[1]])
            if peaks[0] > 100:
                dangling_subreads.append(seq[:peaks[0]])
                qual_dangling_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100:
                dangling_subreads.append(seq[peaks[-1]:])
                qual_dangling_subreads.append(qual[peaks[-1]:])
        else:
            dangling_subreads.append(seq[:peaks[0]])
            qual_dangling_subreads.append(qual[:peaks[0]])
            dangling_subreads.append(seq[peaks[0]:])
            qual_dangling_subreads.append(qual[peaks[0]:])


        # Calculate repeats
        repeats = len(subreads)
        
        # Return subreads data for batched processing instead of calling determine_consensus
        consensus_data = {
            'read': read,
            'subreads': subreads,
            'qual_subreads': qual_subreads,
            'dangling_subreads': dangling_subreads,
            'qual_dangling_subreads': qual_dangling_subreads,
            'repeats': repeats,
            'use': use
        }
    else:
        consensus_data = None
    
    # Collect timing stats (only if variables are defined)
    timing_info = {}
    if 'conk_time' in locals():
        timing_info['conk'] = conk_time
    if 'peaks_time' in locals():
        timing_info['peaks'] = peaks_time  
    
    return consensus_data,read_adapter,subs,peaks,timing_info

def batch_consensus_generation(consensus_batch, abpoa, tmp_dir, args):
    """
    Batch process multiple reads through abpoa for efficiency.
    Returns dict mapping read_name -> consensus_sequence
    """
    if not consensus_batch:
        return {}
    
    # Create temporary files for each read's subreads
    batch_files = []
    read_names = []
    
    for i, (index, consensus_data, adapter, subs) in enumerate(consensus_batch):
        # Show progress every batch_size reads
        if (i + 1) % args.batch_size == 0:
            print(f'\r\tConsensus: {i + 1}/{len(consensus_batch)} subread groups processed...', end='', flush=True)
        read_name = consensus_data['read'][0]
        subreads = consensus_data['subreads']
        qual_subreads = consensus_data['qual_subreads']
        repeats = consensus_data['repeats']
        
        # Skip single subreads (no abpoa needed)
        if repeats == 1:
            continue
            
        # Create FASTQ file for this read's subreads
        subread_file = f'{tmp_dir}batch_read_{i}_{read_name.replace("/", "_")}.fastq'
        with open(subread_file, 'w') as f:
            for j, (subread, qual) in enumerate(zip(subreads, qual_subreads)):
                f.write(f'@{read_name}_subread_{j+1}\n{subread}\n+\n{qual}\n')
        
        batch_files.append(subread_file)
        read_names.append(read_name)
    
    consensus_results = {}
    
    if batch_files:
        # Create list file for abpoa
        list_file = f'{tmp_dir}batch_list.txt'
        with open(list_file, 'w') as f:
            for batch_file in batch_files:
                f.write(f'{batch_file}\n')
        
        # Run abpoa in batch mode
        batch_output = f'{tmp_dir}batch_output.fasta'
        
        # Determine abpoa parameters (use conservative settings)
        insert_lengths = []
        for index, consensus_data, adapter, subs in consensus_batch:
            if consensus_data['repeats'] > 1:
                insert_lengths.extend([len(sub) for sub in consensus_data['subreads']])
        
        if insert_lengths:
            max_length = max(insert_lengths)
            if max_length < 8000:
                abpoa_cmd = f'{abpoa} -M 5 -r 0 -l {list_file} > {batch_output} 2>> abpoa.messages'
            else:
                abpoa_cmd = f'{abpoa} -M 5 -r 0 -S -l {list_file} > {batch_output} 2>> abpoa.messages'
            
            os.system(abpoa_cmd)
            
            # Parse results back to reads
            if os.path.exists(batch_output) and os.path.getsize(batch_output) > 0:
                consensus_seqs = []
                for name, seq, qual in mm.fastx_read(batch_output):
                    consensus_seqs.append(seq)
                
                # Map consensus sequences back to read names
                for i, read_name in enumerate(read_names):
                    if i < len(consensus_seqs):
                        consensus_results[read_name] = consensus_seqs[i]
            
            # Clean up batch files
            try:
                os.remove(list_file)
                os.remove(batch_output)
                for batch_file in batch_files:
                    os.remove(batch_file)
            except:
                pass
    
    # Handle single-subread cases (just use the subread as consensus)
    for index, consensus_data, adapter, subs in consensus_batch:
        read_name = consensus_data['read'][0]
        if consensus_data['repeats'] == 1 and read_name not in consensus_results:
            consensus_results[read_name] = consensus_data['subreads'][0]
    
    return consensus_results

def batch_racon_polishing(racon_batch, racon, tmp_dir, args):
    """
    Batch process multiple consensuses through racon for polishing.
    Returns dict mapping read_name -> polished_consensus_sequence
    """
    if not racon_batch:
        return {}
    
    batch_results = {}
    
    # Group reads into smaller batches for racon processing  
    batch_size = min(args.batch_size, len(racon_batch))  # Use same batch size as abpoa
    
    for batch_start in range(0, len(racon_batch), batch_size):
        batch_end = min(batch_start + batch_size, len(racon_batch))
        current_batch = racon_batch[batch_start:batch_end]
        
        if not current_batch:
            continue
            
        # Create batch files for this sub-batch
        batch_consensus_file = f'{tmp_dir}batch_consensus_{batch_start}.fasta'
        batch_subreads_file = f'{tmp_dir}batch_subreads_{batch_start}.fastq'
        batch_overlaps_file = f'{tmp_dir}batch_overlaps_{batch_start}.paf'
        batch_output_file = f'{tmp_dir}batch_racon_output_{batch_start}.fasta'
        
        # Write all consensuses and subreads for this batch
        with open(batch_consensus_file, 'w') as cons_fh, \
             open(batch_subreads_file, 'w') as sub_fh, \
             open(batch_overlaps_file, 'w') as overlap_fh:
            
            for j, (read_name, consensus_data, abpoa_consensus) in enumerate(current_batch):
                # Show progress every batch_size consensuses  
                global_idx = batch_start + j + 1
                if global_idx % args.batch_size == 0:
                    print(f'\r\tPolishing: {global_idx}/{len(racon_batch)} consensuses processed...', end='', flush=True)
                    
                subreads = consensus_data['subreads']
                qual_subreads = consensus_data['qual_subreads']
                dangling_subreads = consensus_data['dangling_subreads']
                qual_dangling_subreads = consensus_data['qual_dangling_subreads']
                
                # Write consensus sequence
                cons_fh.write(f'>{read_name}\n{abpoa_consensus}\n')
                
                # Create aligner for this consensus
                mm_align = mm.Aligner(seq=abpoa_consensus, preset='map-ont')
                
                # Write subreads and generate overlaps
                for i, (subread, qual) in enumerate(zip(subreads, qual_subreads)):
                    qname = f'{read_name}_{i+1}'
                    sub_fh.write(f'@{qname}\n{subread}\n+\n{qual}\n')
                    
                    # Generate overlap for this subread
                    for hit in mm_align.map(subread):
                        overlap_fh.write(f'{qname}\t{len(subread)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t{read_name}\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}\n')
                
                # Write dangling subreads and generate overlaps
                for j, (subread, qual) in enumerate(zip(dangling_subreads, qual_dangling_subreads)):
                    if j == 0:
                        qname = f'{read_name}_{j}'
                    else:
                        qname = f'{read_name}_{len(subreads) + j}'
                    sub_fh.write(f'@{qname}\n{subread}\n+\n{qual}\n')
                    
                    # Generate overlap for this dangling subread
                    for hit in mm_align.map(subread):
                        overlap_fh.write(f'{qname}\t{len(subread)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t{read_name}\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}\n')
        
        # Run racon on this batch
        racon_msgs_file = f'{tmp_dir}batch_racon_messages_{batch_start}.log'
        
        # Calculate appropriate window size based on insert lengths
        insert_lengths = []
        for read_name, consensus_data, abpoa_consensus in current_batch:
            for subread in consensus_data['subreads']:
                insert_lengths.append(len(subread))
        
        if insert_lengths:
            window = str(min(int(np.median(insert_lengths) * 1.5), 3000))
        else:
            window = "1500"  # default
        
        # Run racon command
        racon_cmd = [racon, batch_subreads_file, batch_overlaps_file, batch_consensus_file, 
                    '-q', '5', '-t', str(args.numThreads), '-w', window]
        
        try:
            with open(batch_output_file, 'w') as output_fh, open(racon_msgs_file, 'w') as msgs_fh:
                subprocess.run(racon_cmd, stdout=output_fh, stderr=msgs_fh, check=True)
        except subprocess.CalledProcessError:
            # If racon fails, fall back to abpoa consensus
            for read_name, consensus_data, abpoa_consensus in current_batch:
                batch_results[read_name] = abpoa_consensus
            continue
        
        # Parse racon results
        if os.path.exists(batch_output_file) and os.path.getsize(batch_output_file) > 0:
            for name, seq, qual in mm.fastx_read(batch_output_file):
                batch_results[name] = seq
        else:
            # Fall back to abpoa consensus if racon output is empty
            for read_name, consensus_data, abpoa_consensus in current_batch:
                batch_results[read_name] = abpoa_consensus
        
        # Clean up batch files
        try:
            for file_path in [batch_consensus_file, batch_subreads_file, batch_overlaps_file, 
                            batch_output_file, racon_msgs_file]:
                if os.path.exists(file_path):
                    os.remove(file_path)
        except:
            pass
    
    return batch_results

def finish_consensus_with_racon(args, consensus_data, racon_consensus, racon, tmp_dir):
    """
    Complete consensus generation with racon polishing results
    """
    read = consensus_data['read']
    subreads = consensus_data['subreads'] 
    repeats = consensus_data['repeats']
    
    name, seq, qual = read[0], read[1], read[2]
    
    # Create subs list
    subs = []
    for i, subread in enumerate(subreads):
        qname = name + '_' + str(i+1)
        subs.append((qname, subread, consensus_data['qual_subreads'][i]))
    
    # Add dangling subreads
    for j, subread in enumerate(consensus_data['dangling_subreads']):
        if j == 0:
            qname = name + '_' + str(j)
        else:
            qname = name + '_' + str(len(subreads) + j)
        subs.append((qname, subread, consensus_data['qual_dangling_subreads'][j]))
    
    if racon_consensus:
        # Calculate quality score
        numbers = []
        for Q in qual:
            numbers.append(ord(Q)-33)
        avg_qual = str(round(np.average(numbers), 1))
        seq_len = len(seq)
        cons_len = len(racon_consensus)
        
        consensus = f'>{name}_{avg_qual}_{seq_len}_{repeats}_{cons_len}\n{racon_consensus}\n'
        return consensus, repeats, subs
    
    return '', repeats, subs

def create_files(adapter,args,outDict,outSubDict,outCountDict):
    outCountDict[adapter]=0
    if not os.path.exists(args.out_path + adapter):
        os.mkdir(args.out_path + adapter)
    outCountDict[adapter]=0
    if args.compress_output:
        writeMode='ab+'
        outDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Consensus.fasta.gz',writeMode)
        outSubDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Subreads.fastq.gz',writeMode)
    else:
        writeMode='a'
        outDict[adapter]=open(args.out_path + adapter +'/R2C2_Consensus.fasta',writeMode)
        outSubDict[adapter]=open(args.out_path + adapter +'/R2C2_Subreads.fastq',writeMode)
    return outDict,outSubDict,outCountDict


def main(args):
    if not args.out_path.endswith('/'):
        args.out_path += '/'
    log_file = open(args.out_path + 'c3poa.log', 'a+')
    
    if not args.nosplint:
        splint_dict = {}
        for splint in mm.fastx_read(args.splint_file, read_comment=False):
            splint_dict[splint[0]] = [splint[1]]
            splint_dict[splint[0]].append(mm.revcomp(splint[1]))

    tmp_dir = args.out_path + 'tmp/'
    reads=args.reads
    fileStart=time.time()

    outDict={}
    outSubDict={}
    outCountDict={}
    previous=set()

    total_reads=0
    short_reads=0
    no_splint_reads=0
    consNumber=0
    
    # Start preprocessing timing (includes minimap2 + task prep + conk/peaks)
    preprocessing_start = time.time()
    
    if not args.nosplint:
        print(f'\tStarting minimap2 alignment...')
        minimap2_start = time.time()
        adapter_dict, adapter_set, no_splint = preprocess(args, tmp_dir, reads)
        minimap2_end = time.time()
        minimap2_total = minimap2_end - minimap2_start
        print(f'\tPreprocessing (minimap2): completed in {minimap2_total:.2f}s')
        for adapter in adapter_set:
            outDict,outSubDict,outCountDict=create_files(adapter,args,outDict,outSubDict,outCountDict)
    else:
        adapter_dict={}
        adapter_set=set()
        adapter_set.add('noSplint')
        for name,seq,q in mm.fastx_read(reads):
            adapter_dict[name]=seq[200:400]
        outDict,outSubDict,outCountDict=create_files('noSplint',args,outDict,outSubDict,outCountDict)
    
    # Prepare tasks for concurrent processing
    tasks = []
    processed_file=open(f'{reads}_processed','w')
    read_count = 0
    for name,seq,q in mm.fastx_read(reads, read_comment=False):
        read_count += 1
        total_reads+=1
        if len(seq) < args.lencutoff:
             short_reads+=1
             processed_file.write(f'{name}\n')
        elif name not in adapter_dict:
             no_splint_reads+=1
             processed_file.write(f'{name}\n')
        else:
            if not args.nosplint:
                read_adapter_info=adapter_dict[name]
                strand = read_adapter_info[1]
                if strand == '-':
                    # use reverse complement of the splint
                    splint = splint_dict[read_adapter_info[0]][1]
                else:
                    splint = splint_dict[read_adapter_info[0]][0]
                tasks.append((name, args, [name,seq,q], splint, read_adapter_info[0], racon, tmp_dir, abpoa))
            else:
                 splint=adapter_dict[name]
                 tasks.append((name, args, [name,seq,q], splint, 'noSplint', racon, tmp_dir, abpoa))

    print(f'\tStarting preprocessing for {len(tasks)} reads...')

    # Process with concurrent.futures for better progress tracking
    results = {}
    completed_count = 0
    
    with ProcessPoolExecutor(max_workers=args.numThreads) as executor:
        # Submit all tasks
        future_to_name = {
            executor.submit(analyze_reads, task[1], task[2], task[3], task[4], task[5], task[6], task[7]): task[0] 
            for task in tasks
        }
        
        # Process completed tasks with progress feedback
        last_progress_time = time.time()
        for future in as_completed(future_to_name):
            current_time = time.time()
            
            name = future_to_name[future]
            results[name] = future.result()
            completed_count += 1
            
            # Show progress every batch_size completions OR every 10 seconds
            if (completed_count % args.batch_size == 0 or 
                completed_count == len(tasks) or
                current_time - last_progress_time > 10):
                last_progress_time = current_time
                print(f'\r\tPreprocessing (conk + peaks): {completed_count}/{len(tasks)} reads processed...', end='', flush=True)
    
    preprocessing_end = time.time()
    preprocessing_total = preprocessing_end - preprocessing_start
    print(f'\r\tPreprocessing (conk + peaks): {len(tasks)}/{len(tasks)} reads completed in {preprocessing_total:.2f}s')
    gc.collect()


    # Phase 1: Collect all results 
    adapters_with_reads=set()
    timing_stats = {'conk': [], 'peaks': [], 'abpoa': [], 'racon': []}
    consensus_batch = []  # List of (index, consensus_data, adapter) for batch processing
    
    for index, result in results.items():
        consensus_data, adapter, subs, peaks, timing_info = result
        
        # Collect timing stats
        for operation, duration in timing_info.items():
            if operation in timing_stats:
                timing_stats[operation].append(duration)
        
        if consensus_data:
            # Add to batch for processing
            consensus_batch.append((index, consensus_data, adapter, subs))
        else:
            # No consensus data (read was filtered out)
            processed_file.write(f'{index}\n')
    
    # Phase 2: Batch process all consensus generation
    print(f'\tStarting abPOA consensus generation for {len(consensus_batch)} subread groups...')
    batch_start_time = time.time()
    
    if consensus_batch:
        consensus_results = batch_consensus_generation(consensus_batch, abpoa, tmp_dir, args)
        
        batch_time = time.time() - batch_start_time
        timing_stats['abpoa'].append(batch_time)
        print(f'\r\tConsensus: {len(consensus_batch)}/{len(consensus_batch)} subread groups completed in {batch_time:.2f}s')
        
        # Phase 3: Racon polishing
        racon_batch = []
        for index, consensus_data, adapter, subs in consensus_batch:
            read_name = consensus_data['read'][0]
            
            if read_name in consensus_results:
                abpoa_consensus = consensus_results[read_name]
                racon_batch.append((read_name, consensus_data, abpoa_consensus))
        
        print(f'\tStarting racon polishing for {len(racon_batch)} consensuses...')
        racon_start_time = time.time()
        
        racon_results = {}
        if racon_batch:
            racon_results = batch_racon_polishing(racon_batch, racon, tmp_dir, args)
            
            racon_time = time.time() - racon_start_time
            timing_stats['racon'].append(racon_time)
            print(f'\r\tPolishing: {len(racon_batch)}/{len(racon_batch)} consensuses completed in {racon_time:.2f}s')
        
        # Phase 4: Write final results
        for index, consensus_data, adapter, subs in consensus_batch:
            read_name = consensus_data['read'][0]
            
            if read_name in consensus_results:
                abpoa_consensus = consensus_results[read_name]
                
                # Get racon polished consensus or fall back to abpoa consensus
                racon_consensus = racon_results.get(read_name, abpoa_consensus)
                
                # Complete consensus generation
                final_consensus, repeats, final_subs = finish_consensus_with_racon(
                    args, consensus_data, racon_consensus, racon, tmp_dir
                )
                
                if final_consensus:
                    if args.compress_output:
                        final_consensus = final_consensus.encode()
                    outDict[adapter].write(final_consensus)
                    outCountDict[adapter] += 1
                    consNumber += 1
                    
                    # Write subreads
                    for subname, subseq, subq in final_subs:
                        entry = f'@{subname}\n{subseq}\n+\n{subq}\n'
                        if args.compress_output:
                            entry = entry.encode()
                        outSubDict[adapter].write(entry)
            
            processed_file.write(f'{index}\n')
    processed_file.close()
    log_file.write(f'Too short reads: {short_reads}'+' ({:.2f}%)'.format((short_reads/total_reads)*100)+'\n')
    log_file.write(f'No splint reads: {no_splint_reads}'+' ({:.2f}%)'.format((no_splint_reads/total_reads)*100)+'\n')
    log_file.write(f'Successful consensus reads: {consNumber}'+' ({:.2f}%)'.format((consNumber/total_reads)*100)+'\n')

    fileEnd=time.time()
    duration=fileEnd-fileStart
    
    
    print(f'\tFinished generating {"{:,}".format(consNumber)} consensus sequences from {"{:,}".format(total_reads)} raw reads ({round((consNumber/total_reads)*100)}%).')
    for adapter in adapter_set:
        outDict[adapter].close()
        outSubDict[adapter].close()
        log_file.write(f'\t{outCountDict[adapter]} consensus reads generated for {adapter}\n')
    log_file.flush()
    log_file.close()


if __name__ == '__main__':
    args = parse_args()
    mp.set_start_method("spawn")
    main(args)

