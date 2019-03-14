# fetch_reads
fetch_reads.py in falcon_unzip is very slow, so i re-write it using C language

fetch_reads.py is a python script file of falcon-unzip, falcon-unzip was included in pb-assembly

setp1:
find file "pb-assembly/lib/python2.7/site-packages/falcon_unzip/tasks/unzip.py",open it using VIM

Comment the line below " # Those outputs are used only by fetch_reads."
#python -m falcon_unzip.mains.fetch_reads --p-ctg={input.p_ctg} --fofn={input.fofn} --ctg-list={output.ctg_list_file}
ulimit -n 4096

Add line
PATH/fetch_reads -p {input.p_ctg} -f {input.fofn} -c {output.ctg_list_file}

Save it

g++ -o fetch_reads fetch_reads.cpp -I/share/public/software/boost1.57.0/include  -L/share/public/software/boost1.57.0/lib/ -lboost_regex
