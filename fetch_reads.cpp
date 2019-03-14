#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <map>
#include <vector>
#include <sys/stat.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string_regex.hpp>


using namespace::std;

void usage(){
	cerr << "\n\nfetch_reads [options]: \n";
	cerr << "  --ctg-list-fn    CTG_LIST_FN   Output. List of ctg_id names. (default: ctg_list)\n";
	cerr << "  --base-dir       BASE_DIR      the base working dir of a falcon assembly (default: ../..\n";
	cerr << "  --p-ctg-fn       P_CTG_FN      Primary contigs of falcon assembly. (FASTA) (default: {base_dir}/2-asm-falcon/p_ctg.fa)\n";
	cerr << "  --fofn           FOFN          path to the file of the list of raw read fasta files (default: ../../input.fofn) \n";
	cerr << "  --ctg-id         CTG_ID        contig identifier in the contig fasta file (default: all)\n";
	cerr << "  --min-ctg-length MIN_CTG_LENTH the minimum length of the contig for the outputs (default: 20000)\n";
	cerr << "Using the read to contig mapping data, to partition the reads (into {ctg_id}_reads.fa and {ctg_id}_ref.fa) grouped by contigs.\n\n\n";
	cerr << "Author: Haikuan Zhang,Mail:zhk_@hotmail.com\n\n\n";
}

// input ">000000F 000011093:B~000043035:E~000016365:B~000194748:E ctg_linear" drop the useless string and only output "000000F"
string GetName(string &line){
	boost::regex equal("\\s+");
	boost::sregex_token_iterator it(line.begin(), line.end(), equal, -1);
	boost::sregex_token_iterator end;
	if(it!=end){
		//cout << *it << endl;
		return *it;
	}else{
		cerr << "ERROR: function GetName..." << endl;
	}
}

std::vector<std::string> SplitStr(std::string str, char delim){
	std::vector<std::string> vs;
	std::string temp;
	std::string new_str = str;
	for(int i=0; i<new_str.length(); i++){
		if(new_str[i]!=delim){
			temp+=new_str[i];
		}else{
			if(temp.compare("") != 0){ vs.push_back(temp);}
			temp = "";
		}
	}
	if(temp.compare("") != 0){ vs.push_back(temp);}
	return vs;
}

// if a file is empty return ture, else false
bool isemptyfile(string filepath){
	ifstream file(filepath.c_str());
	if(!file){
		return true;
	}
	int c=file.get();
	if(file.eof()){
		file.close();
		return true;
	}else{
		file.close();
		return false;
	}
}

int main(int argc, char *argv[]){
	if(argc<=1){
		usage();
		return 1;
	}
	//varable of parameters
	string ctg_list_fn("ctg_list");
	string base_dir("../..");
	string p_ctg_fn(base_dir+"/2-asm-falcon/p_ctg.fa");
	string fofn("../../input.fofn");
	string ctg_id("all");
	unsigned long long min_ctg_length = 20000;
	bool help = false;
	
	const char * shortOpt = "c:b:p:f:g:m:h";
	struct option longOpt[]={
								{"ctg-list-fn", 1, NULL, 'c'},
								{"base-dir", 1, NULL, 'b'},
								{"p-ctg-fn", 1, NULL, 'p'},
								{"fofn", 1, NULL, 'f'},
								{"ctg-id", 1, NULL, 'g'},
								{"min-ctg-length", 1, NULL, 'm'},
								{"help", 0, NULL, 'h'}
							};
	int nextOpt;
	while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1){
		switch (nextOpt){
			case 'c':
				ctg_list_fn = optarg;
				break;
			case 'b':
				base_dir = optarg;
				break;
			case 'p':
				p_ctg_fn = optarg;
				break;
			case 'f':
				fofn = optarg;
				break;
			case 'g':
				ctg_id = optarg;
				break;
			case 'm':
				min_ctg_length = atoi(optarg);
				break;
			case 'h':
				help = true;
				break;
			default:
				usage();
				return 0;
		}
	} 

//Start do analysis	
	string ctg_fa(p_ctg_fn);
	string out_dir(".");
	string read_map_dir(".");
	string rawread_id_file(read_map_dir+"/dump_rawread_ids/rawread_ids");
	string pread_id_file(read_map_dir+"/dump_pread_ids/pread_ids");

	map<string, FILE *> list;

	vector<string> rid_to_oid;
	vector<string> pid_to_fid;
	string line;
	ifstream rtoo(rawread_id_file.c_str());
	while(getline(rtoo, line)){
		rid_to_oid.push_back(line);
	}
	rtoo.close();

	ifstream ptoo(pread_id_file.c_str());
	while(getline(ptoo, line)){
		pid_to_fid.push_back(line);
	}
	ptoo.close();
	

	ifstream ref_fasta(ctg_fa.c_str());
	string chr("");
	string seq("");
	string old_chr("");
	while(getline(ref_fasta, line)){
		if(line == ""){cout << "Blank line!" << endl;continue;}
		if(line[0] == '>'){
			chr = GetName(line);

			if(old_chr != ""){
				if(!ctg_id.compare("all") && !old_chr.compare(ctg_id)){
					old_chr = chr;
					continue;
				}
				if(seq.length()<min_ctg_length){
					old_chr = chr;
					continue;
				}
			
				
				string dir= out_dir+"/"+old_chr.substr(1);
				mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				ofstream ref_out(string(dir+"/ref.fa").c_str());
				FILE *of=fopen(string(dir+"/reads.fa").c_str(), "w");
				if(of==NULL){
					cout << "couldn't open file " << dir << "/reads.fa" << endl;
				}
				list[old_chr.substr(1)] = of;
				ref_out << old_chr << endl;
				ref_out << seq << endl;
				seq="";
				ref_out.close();
			}
			old_chr = chr;
		}else{
			seq = seq+line;
		}
	}
	if((ctg_id == "all" || ctg_id==old_chr) && seq.length()>=min_ctg_length){
		string dir= out_dir+"/"+old_chr.substr(1);
		mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		FILE *of=fopen(string(dir+"/reads.fa").c_str(), "w");
		if(of==NULL){
			cout << "couldn't open file " << dir << "/reads.fa" << endl;
		}
		list[old_chr.substr(1)] = of;
		
		ofstream ref_out(string(dir+"/ref.fa").c_str());
		ref_out << old_chr << endl;
		ref_out << seq << endl;
		ref_out.close();
	}
	ref_fasta.close();
	string u_dir = out_dir+"/unassigned";
	mkdir(u_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	FILE *ofu=fopen(string(u_dir+"/reads.fa").c_str(), "w");
	list["unassigned"] = ofu;

	map<string, string> read_set;
	map<string, long long> ctg_id_hits;
	string map_fn_raw(read_map_dir+"/rawread_to_contigs");
	string map_fn_prd(read_map_dir+"/pread_to_contigs");
	ifstream f_raw(map_fn_raw.c_str());
	ifstream f_prd(map_fn_prd.c_str());
	vector<string> row;
	while(getline(f_raw, line)){
		row = SplitStr(line, ' ');
		string hit_ctg = SplitStr(row[1],'-')[0];
		if(int(atoi(row[3].c_str()))==0){
			string o_id = rid_to_oid[int(atoi(row[0].c_str()))];
			read_set[o_id] = hit_ctg;
			ctg_id_hits[hit_ctg]++;
		}
	}
	map<string, string>::iterator itt;
	while(getline(f_prd, line)){
		row = SplitStr(line, ' ');
		string hit_ctg = SplitStr(row[1], '-')[0];
		itt = read_set.find(hit_ctg);
		if(itt == read_set.end() && int(atoi(row[3].c_str()))==0){
			row = SplitStr(pid_to_fid[atoi(row[0].c_str())],'/');
			int index = atoi(row[1].c_str())/10;
			string o_id = rid_to_oid[index];
			read_set[o_id] = hit_ctg;
			ctg_id_hits[hit_ctg]++;
		}
	}
	f_raw.close();
	f_prd.close();
	
	ifstream rfn(fofn.c_str());
	string line2;
	while(getline(rfn, line)){
		ifstream read_fa_file(line.c_str());
		string reads_id("");
		string reads_seq("");
		string old_reads_id("");
		while(getline(read_fa_file, line2)){
			if(line2 == ""){cout << "Blank line!" << endl;continue;}
			if(line2[0] == '>'){
				reads_id = GetName(line2);
				if(old_reads_id != ""){
			//		cout << old_reads_id << endl;
			//		cout << reads_seq << endl;
			//		cout << old_reads_id.substr(1) << endl;
					itt = read_set.find(old_reads_id.substr(1));
					if(itt != read_set.end() && list.find(read_set[old_reads_id.substr(1)]) != list.end()){
						fprintf(list[read_set[old_reads_id.substr(1)]], "%s\n%s\n", old_reads_id.c_str(), reads_seq.c_str());
					}else{
						fprintf(list["unassigned"], "%s\n%s\n", old_reads_id.c_str(), reads_seq.c_str());
					}
					reads_seq = "";
				}
				old_reads_id = reads_id;
			}else{
				reads_seq+=line2;
			}
		}
		itt = read_set.find(old_reads_id.substr(1));
		if(itt != read_set.end() && list.find(read_set[old_reads_id.substr(1)]) != list.end()){
			fprintf(list[read_set[old_reads_id.substr(1)]], "%s\n%s\n", old_reads_id.c_str(), reads_seq.c_str());
		}else{
			fprintf(list["unassigned"], "%s\n%s\n", old_reads_id.c_str(), reads_seq.c_str());
		}
	}
	map<string, FILE *>::iterator it_h = list.begin();
	while(it_h != list.end()){
		fclose(it_h->second);
		it_h++;
	}

	//Log 
	map<string, long long>::iterator it_hits=ctg_id_hits.begin();
	ofstream log(string(out_dir+"/"+ctg_list_fn).c_str());
	while(it_hits != ctg_id_hits.end()){
		if(it_hits->first[it_hits->first.size()-1] != 'F' && it_hits->first[it_hits->first.size()-1] != 'R'){
			cout << "Skipping contig: {}. Reason: ignoring small circular contigs." << it_hits->first << endl;
		}else if(it_hits->second<5){
			cout << "Skipping contig: {}. Reason: not enough reads (ctg_id_hits.get(ctg_id, 0) = {} < 5)." << it_hits->first << endl;
		}else if(isemptyfile(out_dir+"/"+it_hits->first+"/ref.fa") || isemptyfile(out_dir+"/"+it_hits->first+"/reads.fa")){
			cout << "Skipping contig: {}. Reason: non-existent or empty reads/ref file ({})." << it_hits->first   << endl;
		}else{
			log << it_hits->first << endl;
		}
		it_hits++;
	}
	log.close();

/*
	ifstream input("TTT.txt");
	string file;
	string value;
	map<string, FILE *> list;
	map<string, FILE *>::iterator it;
	while(getline(input, line)){
		istringstream istring(line);
		istring >> file;
		istring >> value;
		it=list.find(file);
		if(it != list.end()){
			fprintf(list[file], "%s\t%d\n", value.c_str(), 10);
		}else{
			FILE *of=fopen(file.c_str(), "w");
			list[file] = of;
			fprintf(list[file], "%s\t%d\n", value.c_str(), 10);
		}
	}
	it = list.begin();
	while(it!=list.end()){
		fclose(it->second);
		it++;
	}
*/
	return 0;
}
