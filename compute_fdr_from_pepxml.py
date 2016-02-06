#!/bin/python
#
# compute_fdr_from_pepxml.py - parse any number of pepxml files matching a glob
# pattern and compute fdr and q-value, then write to an output file.
#
# Uses linear combination of comet search engine scores by Philip Wilmarth,
# based on the original method published by Keller, et al (2002).

# standard libary imports...
import sys
import math
import glob
import argparse
import os.path
import time
import multiprocessing
import os
import gc
import zipfile
import re
from Queue import Queue, Empty, Full
from threading import Thread
from itertools import chain
from collections import defaultdict

# local custom imports...
from XMLStackParserForExpat import XMLStackParser, xml_stack_start_handler, xml_stack_end_handler


class ProgressCounter(object):
    def __init__(self, report_interval, report_template):
        self.report_interval = report_interval
        self.report_template = report_template
        self.value = 0

    def __call__(self, thing):
        return self.count(thing)

    def count(self, thing):
        self.value += 1
        if (self.value % self.report_interval) == 0:
            self.report()
        return thing

    def report(self, report_template = None):
        if report_template:
            sys.stderr.write(report_template % self.value)
        else:
            sys.stderr.write(self.report_template % self.value)
        sys.stderr.flush()
        
    def message(self, msg):
        sys.stderr.write(msg)
        sys.stderr.flush()


def compute_new_dcn_from_psm(psm, hit_index):
    # deltaCN is the difference between the current and next best
    # xcorr value, divided by the current xcorr value.
    # deltaCN is computed by coment and included in the pepxml output.
    #
    # newDeltaCN is the difference between the current and average
    # xcorr value, divided by the current xcorr value, where the
    # average xcorr value is defined as the average xcorr of the
    # 4th rank hit (index = 3) all the way up to the top hit (e.g. the 12th hit)
    xcorr_hit = float(psm['search_hits'][hit_index]['search_scores']['xcorr'])
    bottom_xcorrs = [float(hit['search_scores']['xcorr']) for hit in psm['search_hits'][3:]]
    sum_bottom_xcorrs = sum(bottom_xcorrs)
    if sum_bottom_xcorrs and xcorr_hit > 0.0:
        # use the newdcn...
        return (xcorr_hit - sum_bottom_xcorrs/len(bottom_xcorrs))/xcorr_hit
    else:
        # just return the old dcn...
        return float(psm['search_hits'][hit_index]['search_scores']['deltacn'])


def phillips_hit_score(psm, hit_index):
    assumed_charge = int(psm['spectrum_query']['assumed_charge'])
    hit = psm['search_hits'][hit_index]
    protein = hit['protein']
    peptide = hit['peptide']
    peptide_length = len(peptide) if len(peptide) > 2 else 2
    scores = hit['search_scores']
    mass_diff = float(hit['massdiff'])
    xcorr = float(scores['xcorr'])
    if xcorr < 0.01:
        xcorr = 0.01
    sprank = int(scores['sprank'])
    # deltacn = float(scores['deltacn'])    
    deltacn = compute_new_dcn_from_psm(psm, hit_index)
                
    if assumed_charge == 1:
        peptide_length = peptide_length if peptide_length < 100 else 100
        num_ions = 2*peptide_length
        xcorr_prime = math.log(xcorr)/math.log(num_ions)
        phillips_score = (10.0*xcorr_prime +
                7.0*deltacn +
                (-0.25 * math.log(sprank)) +
                (-0.45 * abs(mass_diff)) +
                (-0.05))
    elif assumed_charge == 2:
        peptide_length = peptide_length if peptide_length < 15 else 15
        num_ions = 2*peptide_length
        xcorr_prime = math.log(xcorr)/math.log(num_ions)
        phillips_score = (8.0*xcorr_prime +
                7.5*deltacn +
                (-0.20 * math.log(sprank)) +
                (-0.55 * abs(mass_diff)) +
                (-0.80))
    elif assumed_charge == 3:
        peptide_length = peptide_length if peptide_length < 25 else 25
        num_ions = 4*peptide_length
        xcorr_prime = math.log(xcorr)/math.log(num_ions)
        phillips_score = (10.0*xcorr_prime +
                8.0*deltacn +
                (-0.20 * math.log(sprank)) +
                (-0.40 * abs(mass_diff)) +
                (-1.5))
    elif assumed_charge == 4:
        peptide_length = peptide_length if peptide_length < 25 else 25
        num_ions = 4*peptide_length
        xcorr_prime = math.log(xcorr)/math.log(num_ions)
        phillips_score = (10.0*xcorr_prime +
                8.0*deltacn +
                (-0.20 * math.log(sprank)) +
                (-0.40 * abs(mass_diff)) +
                (-1.5))
    else:
        phillips_score = -10.0
    return phillips_score


def is_decoy_protein(protein, decoy_string):
    return False if protein.find(decoy_string) == -1 else True


def parse_pepxml_into_psms(fileobj, label, run_identifier):
    # queue of psms shared by producer and consumer threads...
    psm_queue = Queue()
    
    def producer():
        current_psm = [None]
        with XMLStackParser(fileobj, current_psm) as parser:
            @xml_stack_start_handler('spectrum_query')
            def spectrum_query_start_handler(helper, name, attrs, captures, texts, current_psm):
                psm = {}
                psm['label'] = label
                psm['run_identifier'] = run_identifier
                psm['spectrum_query'] = attrs
                psm['search_hits'] = []
                current_psm[0] = psm
                
            @xml_stack_start_handler('spectrum_query/search_result/search_hit')
            def search_hit_start_handler(helper, name, attrs, captures, texts, current_psm):
                psm = current_psm[0]
                psm['search_hits'].append(attrs)
                current_hit = psm['search_hits'][-1]
                current_hit['search_scores'] = {}
                current_hit['modification_string'] = None
                current_hit['modifications'] = set()
    
            @xml_stack_start_handler('spectrum_query/search_result/search_hit/search_score')
            def search_score_start_handler(helper, name, attrs, captures, texts, current_psm):
                psm = current_psm[0]
                psm['search_hits'][-1]['search_scores'][attrs.name] = attrs.value

            @xml_stack_end_handler('spectrum_query')
            def spectrum_query_end_handler(helper, name, attrs, captures, texts, current_psm):
                psm = current_psm[0]
                psm_queue.put(psm)
                current_psm[0] = None

            @xml_stack_start_handler('spectrum_query/search_result/search_hit/modification_info/mod_aminoacid_mass')
            def search_hit_modification_info_mod_aminoacid_mass_start_handler(helper, name, attrs, captures, texts, current_psm):
                psm = current_psm[0]
                current_hit = psm['search_hits'][-1]
                position = int(attrs.position) - 1
                amino_acid = current_hit['peptide']
                modification = amino_acid[position] + attrs.mass
                current_hit['modifications'].add(modification)                
                current_hit['modification_string'] = '+'.join(psm['search_hits'][-1]['modifications'])
                
                
            # no text to track so turn off text accumulation
            # to save lots of memory...
            parser.accumulate_text(False)
            # now that all the hqndlers are declared,
            # parse the file...
            parser.parse()
            
        # done reading the file...
        fileobj.close()
        # block thread here until all queued items have
        # been processed...
        psm_queue.join()

    # fire off the producer thread as a daemon...
    producer_thread = Thread(target=producer)
    producer_thread.daemon = True
    producer_thread.start()    

    # consumer thread runs as generator,
    # keep running as long as thread is alive,
    # but there has to be a timeout on the blocking
    # get because of the possibility of a synchronization
    # deadlock where all items are processed but the
    # producer thread has not yet had a chance to exit...
    while producer_thread.is_alive():
        try:
            psm = psm_queue.get(block=True, timeout = 1.0)
            psm_queue.task_done()
            yield psm
        except Empty:
            continue
        

def extract_short_filename(filename):
    return os.path.split(filename)[1] 


def extract_filename_extension(filename):
    return filename.split('.')[-1].lower()


def extract_run_identifier(filename, patterns):
    SHORTCUT_PATTERNS = { 'PSR767': "(.*_\d[AB])-?(?:\d)(_.*) \g<1>\g<2>",
                         'FILENAME': "(.*) \g<1>", }
    COMPILED_PATTERN_CACHE = {}

    # allow shortcuts...
    if patterns in SHORTCUT_PATTERNS:
        patterns = SHORTCUT_PATTERNS[patterns]

    # pre-compile pattern since the same
    # pattern is used over and over for the whole run...
    if patterns not in COMPILED_PATTERN_CACHE:
        p0, p1 = patterns.split(' ')
        COMPILED_PATTERN_CACHE[patterns] = (re.compile(p0), p1)
            
    # extract the run identifier from the short file name...
    shortfilename = extract_short_filename(filename)
    r0, p1 = COMPILED_PATTERN_CACHE[patterns]
    run_identifier = r0.sub(p1, shortfilename)
    return run_identifier


def parse_file_into_psms(filename, opts):    
    def is_pep_xml_filename(fname):
        extension = extract_filename_extension(fname)
        if extension == 'xml' or extension == 'pepxml':
            return True
        else:
            return False
    
    def echo_zip_filename(fname):
        shortfilename = extract_short_filename(fname) 
        sys.stderr.write("Unzipping %s..." % shortfilename)
        return fname
    
    
    # open file handles here, but file handles are closed
    # when the XML parsing is completed in parse_pepxml_into_psms()...
    if is_pep_xml_filename(filename):
        fileobj = open(filename, 'rU')
        psms = parse_pepxml_into_psms(fileobj,
                                      extract_short_filename(filename),
                                      extract_run_identifier(filename, opts.run_identifier))
        return psms
    elif extract_filename_extension(filename) == 'zip':
        zf = zipfile.ZipFile(filename, 'r')
        pepxml_fnames = [zfname for zfname in zf.namelist() if is_pep_xml_filename(zfname)]
        psms = chain(*[parse_pepxml_into_psms(zf.open(echo_zip_filename(zfname), 'rU'),
                                              extract_short_filename(zfname),
                                              extract_run_identifier(zfname, opts.run_identifier)) for zfname in pepxml_fnames])        
        zf.close()
        return psms
    else:
        raise ValueError("Unknown file format for '%s'." % filename)


def extract_top_hit(psm):    
    if psm.get('spectrum_query') and psm.get('search_hits') and len(psm['search_hits']) > 0:
        top_score = None
        top_index = None
        for index, hit in enumerate(psm['search_hits']):
            hit_score = phillips_hit_score(psm, index)
            if top_index is None or hit_score > top_score:
                top_score = hit_score
                top_index = index
        top_hit = psm['search_hits'][top_index]
        protein = top_hit['protein']
        peptide = top_hit['peptide']
        filename = psm['label']
        run_identifier = psm['run_identifier']
        start_scan = psm['spectrum_query']['start_scan']
        end_scan = psm['spectrum_query']['end_scan']
        assumed_charge = int(psm['spectrum_query']['assumed_charge'])
        num_tol_term = int(top_hit['num_tol_term'])
        modification_string = top_hit['modification_string']
        return (protein, peptide, filename, run_identifier, start_scan, end_scan, assumed_charge, num_tol_term, modification_string, top_score)
    else:
        return None


def write_output_as_tsv(q_values, outfile):
    headers = ['protein', 'peptide', 'filename', 'run_identifier', 'start_scan', 'end_scan', 'assumed_charge', 'num_tol_term', 'modification_string', 'score', 'fdr', 'qvalue']
    outfile.write("%s\n" % '\t'.join(headers))
    for q in q_values:
        outfile.write("%s\n" % '\t'.join([str(v).replace('\t', ' ') for v in q]))
        
        
def write_output_file(q_values, opts):
    outfile = open(opts.out, 'wt') if opts.out != '-' else sys.stdout
    if opts.out_format == 'tsv':
        write_output_as_tsv(q_values, outfile)
    else:
        raise ValueError("The string '%s' is not a valid output format." % opts.out_format)
    if outfile != sys.stdout:
        outfile.close()


def read_pepxml_file(filename, opts):
    shortfilename = extract_short_filename(filename)
    counter = ProgressCounter(1000, "%d...")
    sys.stderr.write("Reading %s..." % shortfilename)
    top_hits_and_scores = [hit for hit in (extract_top_hit(counter(psm)) for psm in parse_file_into_psms(filename, opts)) if hit]
    counter.report("%d.\n")
    sys.stderr.write("OK.\n")
    return top_hits_and_scores


def read_pepxml_file_into_queue(filename, queue, opts):
    shortfilename = extract_short_filename(filename)
    counter = ProgressCounter(1000, "%%d(%d)..." % os.getpid())
    sys.stderr.write("Reading %s..." % shortfilename)    
    for hit in (extract_top_hit(counter(psm)) for psm in parse_file_into_psms(filename, opts)):
        if hit:
            retry = True
            while retry:
                try:
                    queue.put(hit, block = True, timeout = 1.0)
                    retry = False
                except Full:
                    continue                
    counter.report("%d.\n")
    sys.stderr.write("OK.\n")
    # block thread here until all queued items
    # have been processed...
    queue.join()    
    return True


def read_pepxml_files_multiprocessing(filenames, max_parallel, opts):    
    def create_reader_process(filename, queue):
        reader_process = multiprocessing.Process(target=read_pepxml_file_into_queue, args=(filename, queue, opts))
        reader_process.daemon = True
        return reader_process

    def update_reader_processes(waiting_reader_processes, running_reader_processes):
        for p1 in running_reader_processes:
            if not p1.is_alive():
                running_reader_processes.remove(p1)
                if waiting_reader_processes:
                    p2 = waiting_reader_processes.pop()
                    running_reader_processes.append(p2)
                    p2.start()
            
    # blocking queue to give the reader time to catch up...
    queue = multiprocessing.JoinableQueue(100)
    
    # need to launch a reader process for each file,
    # but no more than one for each CPU core at a time...
    max_running_readers = max_parallel if max_parallel > 0 else multiprocessing.cpu_count()
    waiting_reader_processes = [create_reader_process(filename, queue) for filename in filenames]
    running_reader_processes = []
    for n in range(max_running_readers):
        if waiting_reader_processes:
            p = waiting_reader_processes.pop()
            running_reader_processes.append(p)
            p.start()
            
    # read hits sent by readers via the queue...
    hits = []
    while True:
        try:
            hit = queue.get(block=True, timeout = 5.0)
            hits.append(hit)
            queue.task_done()
            if len(hits) % 1000 == 0:
                # periodically start any waiting readers...
                update_reader_processes(waiting_reader_processes, running_reader_processes)
                # inform user about progress...
                sys.stderr.write("%d(main)..." % len(hits))
        except Empty:
            # nothing to read in queue, start any waiting readers...
            update_reader_processes(waiting_reader_processes, running_reader_processes)
            # if no waiting readers and no running readers
            # then exit loop...
            if not running_reader_processes and not waiting_reader_processes:
                break    
    return hits


def compute_fdr_and_write_output(opts):
    # extract and sort top hits according to highest score is best...
    # read all files that match the glob...
    sys.stderr.write("Reading pepxml files...\n")
    if opts.parallel == 0:
        all_hits_and_scores = list(chain(*(read_pepxml_file(filename, opts) for filename in glob.glob(opts.pepxml))))
    else:
        all_hits_and_scores = read_pepxml_files_multiprocessing(glob.glob(opts.pepxml), opts.max_parallel, opts)
    sys.stderr.write("Finished reading pepxml files. Found %d top hits.\n" % len(all_hits_and_scores))
    
    # separate out FDR computation by peptide class
    # peptide class is defined as (assumed_charge, num_tol_term, num_modifications).
    peptide_classes = defaultdict(lambda : [])
    for h in all_hits_and_scores:
        # peptide hit data structure is:
        # (protein, peptide, filename, run_identifier, start_scan, end_scan, assumed_charge, num_tol_term, modification_string, score)
        protein, peptide, filename, run_identifier, start_scan, end_scan, assumed_charge, num_tol_term, modification_string, score = h
        pc = (assumed_charge, num_tol_term, modification_string)
        peptide_classes[pc].append(h)
        
    # process each peptide class in turn...
    all_q_values = []
    for pc in peptide_classes:        
        # compute fdr values...
        sys.stderr.write("Computing FDR for class %s..." % str(pc))
        sorted_hits_and_scores = sorted(peptide_classes[pc], key=lambda h:h[-1], reverse=True)
        forward_count = 0
        decoy_count = 0
        fdrs = []
        for h in sorted_hits_and_scores:
            if is_decoy_protein(h[0], opts.decoy_string):
                decoy_count += 1
            else:
                forward_count += 1
                fdr = float(decoy_count)/forward_count
                fdrs.append(h + (fdr,))
        sys.stderr.write("OK.\n")
        # compute q-values...
        sys.stderr.write("Computing q-values for class %s..." % str(pc))
        n_fdrs = len(fdrs)
        reversed_q_values = []
        min_q_value = 1.0
        for h in reversed(fdrs):
            min_q_value = min(min_q_value, h[-1])
            reversed_q_values.append(h + (min_q_value,))
        q_values = list(reversed(reversed_q_values))    
        sys.stderr.write("OK.\n")
        # append to full list...
        all_q_values.extend(q_values)
        
    # write results to output...        
    sys.stderr.write("Writing output file...")
    write_output_file(all_q_values, opts)
    sys.stderr.write("OK.\n")
    
    
if __name__ == '__main__':
    # define and parse command line...
    parser = argparse.ArgumentParser(description='Command line for %s.' % sys.argv[0])
    parser.add_argument('-pepxml', type=str, help='Input path or glob of pepXML files or zipped pepXML files to process.', required = True)    
    parser.add_argument('-decoy_string', type=str, help='Character string designating decoy hits. (default = "REV_").', default = "REV_")    
    parser.add_argument('-out', type=str, help='Output file path, or - for stdout. (default = "-").', default = "-")    
    parser.add_argument('-out_format', choices=['tsv',], help='Output file format. (default = "tsv").', default = "tsv")    
    parser.add_argument('-parallel', type=int, help='Flag to turn on reading multiple files in parallel (1 = on, 0 = off, default = 0).', default = 0)    
    parser.add_argument('-max_parallel', type=int, help='Maximum number of simultaneous reader processes. (0 = number cpus, default = 0).', default = 0)    
    parser.add_argument('-run_identifier', type=str, help='Two whitespace-separated regex patterns that extract and substitute run identifiers from input filenames, also one of these shortcuts: {FILENAME, PSR767}.', required = True)    
    opts = parser.parse_args()
    
    # keep track of elapsed time...
    start_time = time.time()
    
    # process files...
    compute_fdr_and_write_output(opts)

    # report elapsed time...
    elapsed_time = time.time() - start_time
    sys.stderr.write("Elapsed time = %0.2f seconds.\n" % elapsed_time)




