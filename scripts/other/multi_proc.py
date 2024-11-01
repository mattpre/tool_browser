from __future__ import print_function
import os,subprocess,sys
import multiprocessing as mp
from queue import Queue
from threading import Thread
import threading
import shutil

def run(f,ws):
    print('Running %s\n' % f)
    ZSoil = r'c:\Program Files\ZSoil\ZSoil 2019 v19.03 x64\Z_Soil.exe '
    subprocess.check_call('%s %s\\%s /E'%(ZSoil,ws,f))

def worker(queue):
    """Process files from the queue."""
    for args in iter(queue.get, None):
        try:
            run(*args)
        except Exception as e: # catch exceptions to avoid exiting the
                               # thread prematurely
            print('%r failed: %s' % (args, e,), file=sys.stderr)

def main():
    # populate files
    ws = r'F:\M1234_Matthieu\coupe_2'
    q = Queue()
    for ks in [1,3,4,0,2]:#range(3):
        for kc,c in enumerate(['Gmin','Gmoy','Gmax']):
            q.put_nowait(('M1234_Coupe2_stab_v1_%s_SET%i.inp'%(c,ks+1),ws))

    # start threads
    threads = [Thread(target=worker, args=(q,)) for _ in range(5)]#mp.cpu_count()-1)]
    for t in threads:
        t.daemon = True # threads die if the program dies
        t.start()
    for _ in threads: q.put_nowait(None) # signal no more files
    for t in threads: t.join() # wait for completion

if __name__ == '__main__':
    mp.freeze_support() # optional if the program is not frozen
    main()
