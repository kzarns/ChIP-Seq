//Created by Adam Burkholder, National Institute of Environmental Health Sciences, 2010-2015

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <sstream>
#include <tr1/unordered_map>
#include <errno.h>
#include <set>
#include <getopt.h>

#ifndef SINGLE
#include <pthread.h>
#endif

using namespace std;
using tr1::unordered_map;

class chr_entry {									//used to store intervals in DB file
	public:
		vector<string> dbdesc1;
		vector<string> dbdesc2;
		vector<long> dbphysical_start;
		vector<long> dbphysical_end;
		vector<string> dbwhole;
};

class chr_entry_s : public chr_entry {				//used for strand-aware matching
	public:
		vector<string> dbstrand;
};

struct option long_options[]={
		{"help",0,NULL,'h'},
		{"strands",1,NULL,'s'},
		{"threads",1,NULL,'t'}
};

struct ptr_entry {									//stores pointers to a single DB file entry
	string *desc1;
	string *desc2;
	string *chr;
	long *physical_start;
	long *physical_end;
	string *strand;
	string *whole;
};

struct totals_info {								//stores most hit locations and counts
	double total;
	double most;
	vector<long> most_loc;
	string chr;
};

ofstream outfile;
ofstream outfile2;
unordered_map<string,chr_entry> db;										//used for strand-independent matching
unordered_map<string,unordered_map<string,chr_entry> > db_split;		//used for same- or opposite-strand only matching
unordered_map<string,chr_entry_s> db_s;									//used for flagged or separate-output sense/antisense matching
map<string,totals_info> table;											//stores most hit info for single output, sense most hit for split
map<string,totals_info> table2;											//stores most hit info for antisense output only
set<string> strand_list;
unordered_map<string,string> strand_map;
int s_val;

#ifndef SINGLE
pthread_mutex_t outlock=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t linelock=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t tablelock=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t tablelock2=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t numlock=PTHREAD_MUTEX_INITIALIZER;
#endif

#ifndef SINGLE
void usage(void) {
	cout << "Usage: cppmatch [options] [DB File Name] [Query File Name] [Output File Name]\n";
	cout << "Available Options:\n";
	cout << "  --help                      produce this help message\n";
	cout << "  -t [ --threads ] arg (=1)   specify number of threads to use\n";
	cout << "  -s [ --strands ] arg (=i)   specify handling of strand identifiers:\n";
	cout << "                                 i    ignore strand identifiers\n";
	cout << "                                 s    report only same-strand hits\n";
	cout << "                                 o    report only opposite-strand hits\n";
	cout << "                                 bf   flag each hit as sense or antisense\n";
	cout << "                                 bs   write sense and antisense hits to separate\n";
	cout << "                                      files\n";
	return;
}
#else
void usage(void) {
	cout << "Usage: cppmatch [options] [DB File Name] [Query File Name] [Output File Name]\n";
	cout << "Available Options:\n";
	cout << "  --help                      produce this help message\n";
	cout << "  -s [ --strands ] arg (=i)   specify handling of strand identifiers:\n";
	cout << "                                 i    ignore strand identifiers\n";
	cout << "                                 s    report only same-strand hits\n";
	cout << "                                 o    report only opposite-strand hits\n";
	cout << "                                 bf   flag each hit as sense or antisense\n";
	cout << "                                 bs   write sense and antisense hits to separate\n";
	cout << "                                      files\n";
	return;
}
#endif

void cache_i(ifstream &infile) {																					//store DB file entries for strand-independent matching
	string line;
	string desc1;
	string desc2;
	string chr;
	long physical_start;
	long physical_end;
	getline(infile,line);
	while(!infile.eof()) {
		istringstream temp1(line);
		temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end;
		if(temp1.fail()) {
			cout << "DB File contains bad line, skipping: " << line << endl;
		}
		else {
			ostringstream temp2;																					//separate entires by chromosome
			temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end;
			db[chr].dbdesc1.push_back(desc1);
			db[chr].dbdesc2.push_back(desc2);
			db[chr].dbphysical_start.push_back(physical_start);
			db[chr].dbphysical_end.push_back(physical_end);
			db[chr].dbwhole.push_back(temp2.str());
			table[desc1].total=0;
			table[desc1].chr=chr;
			table[desc1].most=0;
			table[desc1].most_loc.push_back(0);
		}
		getline(infile,line);
	}
	infile.close();
	return;
}

void cache_split(ifstream &infile) {																							//store DB file entries for same-strand or opposite-strand only matching
	string line;
	string desc1;
	string desc2;
	string chr;
	long physical_start;
	long physical_end;
	string strand;
	getline(infile,line);
	while(!infile.eof()) {
		istringstream temp1(line);
		temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end >> strand;
		if(temp1.fail()) {
			cout << "DB File contains bad line, skipping: " << line << endl;
		}
		else {
			ostringstream temp2;
			temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end << '\t' << strand;
			strand_list.insert(strand);
			db_split[strand][chr].dbdesc1.push_back(desc1);																		//separate entries by strand and chromosome
			db_split[strand][chr].dbdesc2.push_back(desc2);
			db_split[strand][chr].dbphysical_start.push_back(physical_start);
			db_split[strand][chr].dbphysical_end.push_back(physical_end);
			db_split[strand][chr].dbwhole.push_back(temp2.str());
			table[desc1].total=0;
			table[desc1].chr=chr;
			table[desc1].most=0;
			table[desc1].most_loc.push_back(0);
		}
		getline(infile,line);
	}
	infile.close();
	return;
}

void cache_s(ifstream &infile) {																								//store DB file entries for flagged or split output strand-aware matching
	string line;
	string desc1;
	string desc2;
	string chr;
	long physical_start;
	long physical_end;
	string strand;
	getline(infile,line);
	while(!infile.eof()) {
		istringstream temp1(line);
		temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end >> strand;
		if(temp1.fail()) {
			cout << "DB File contains bad line, skipping: " << line << endl;
		}
		else {
			ostringstream temp2;
			temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end << '\t' << strand;
			strand_list.insert(strand);
			db_s[chr].dbdesc1.push_back(desc1);
			db_s[chr].dbdesc2.push_back(desc2);
			db_s[chr].dbphysical_start.push_back(physical_start);
			db_s[chr].dbphysical_end.push_back(physical_end);
			db_s[chr].dbstrand.push_back(strand);
			db_s[chr].dbwhole.push_back(temp2.str());
			table[desc1].total=0;
			table[desc1].chr=chr;
			table[desc1].most=0;
			table[desc1].most_loc.push_back(0);
			table2[desc1].total=0;
			table2[desc1].chr=chr;
			table2[desc1].most=0;
			table2[desc1].most_loc.push_back(0);
		}
		getline(infile,line);
	}
	infile.close();
	return;
}

void query_i(string &line) {																				//perform strand-independent intersection
	bool qstartcmp;
	bool qendcmp;
	string desc1;
	double desc2;
	string chr;
	string whole;
	long physical_start;
	long physical_end;
	ptr_entry arrays;
	ostringstream acc;
	istringstream temp1(line);
	temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end;
	if(temp1.fail()) {
		cout << "Query File contains bad line, skipping: " << line << endl;
	}
	else {
		ostringstream temp2;
		temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end;
		whole=temp2.str();
		if(db.find(chr)!=db.end()) {
			size_t max=db[chr].dbwhole.size();																				//store pointers to first DB entry for current query's chromosome to reduce lookups
			arrays.desc1=&db[chr].dbdesc1[0];
			arrays.desc2=&db[chr].dbdesc2[0];
			arrays.physical_start=&db[chr].dbphysical_start[0];
			arrays.physical_end=&db[chr].dbphysical_end[0];
			arrays.whole=&db[chr].dbwhole[0];
			for(size_t i=0;i<max;i++) {																						//for each entry in query file, check overlap with each entry in DB file
				if(arrays.physical_end[i]<physical_start || arrays.physical_start[i]>physical_end) continue;
				qstartcmp=(physical_start>=arrays.physical_start[i]);														//check for type of overlap and flag
				qendcmp=(physical_end<=arrays.physical_end[i]);
				if(qstartcmp && qendcmp) acc << arrays.whole[i] << '\t' << whole  << "\tB" << endl;
				else if(qstartcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tS" << endl;
				else if(qendcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tE" << endl;
				else acc << arrays.whole[i] << '\t' << whole  << "\tC" << endl;
#ifndef SINGLE
				pthread_mutex_lock(&tablelock);
				totals_info *info=&table[arrays.desc1[i]];																	//check if count in current intersection is greater than or equal to current max
				info->total+=desc2;
				if(desc2>info->most) {
					info->most=desc2;
					info->most_loc.clear();
					info->most_loc.push_back(physical_start);
				}
				else if(desc2==info->most) {
					info->most_loc.push_back(physical_start);
				}
				pthread_mutex_unlock(&tablelock);
			}
			pthread_mutex_lock(&outlock);
			outfile << acc.str();																							//write intersection to output
			pthread_mutex_unlock(&outlock);
#else
				totals_info *info=&table[arrays.desc1[i]];
				info->total+=desc2;
				if(desc2>info->most) {
					info->most=desc2;
					info->most_loc.clear();
					info->most_loc.push_back(physical_start);
				}
				else if(desc2==info->most) {
					info->most_loc.push_back(physical_start);
				}
			}
			outfile << acc.str();
#endif
		}
	}
	return;
}

void query_split(string &line) {																							//perform intersection for same- or opposite-strand only matching
	string desc1;
	double desc2;
	string chr;
	long physical_start;
	long physical_end;
	string strand;
	string whole;
	bool qstartcmp;
	bool qendcmp;
	ptr_entry arrays;
	ostringstream acc;
	istringstream temp1(line);
	temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end >> strand;
	if(temp1.fail()) {
		cout << "Query File contains bad line, skipping: " << line << endl;
	}
	else {
		ostringstream temp2;
		temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end << '\t' << strand;
		whole=temp2.str();
		string str=strand_map[strand];																							//get strand to be matched
		if(str=="") {
			str=strand_map["dummy"];																							//use "dummy" if strand of current query entry was not found in the DB file
		}
		if(db_split.find(str)!=db_split.end()) {
			if(db_split[str].find(chr)!=db_split[str].end()) {
				size_t max=db_split[str][chr].dbwhole.size();
				arrays.desc1=&db_split[str][chr].dbdesc1[0];
				arrays.desc2=&db_split[str][chr].dbdesc2[0];
				arrays.physical_start=&db_split[str][chr].dbphysical_start[0];
				arrays.physical_end=&db_split[str][chr].dbphysical_end[0];
				arrays.whole=&db_split[str][chr].dbwhole[0];
				for(size_t i=0;i<max;i++) {
					if(arrays.physical_end[i]<physical_start || arrays.physical_start[i]>physical_end) continue;
					qstartcmp=(physical_start>=arrays.physical_start[i]);
					qendcmp=(physical_end<=arrays.physical_end[i]);
					if(qstartcmp && qendcmp) acc << arrays.whole[i] << '\t' << whole  << "\tB" << endl;
					else if(qstartcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tS" << endl;
					else if(qendcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tE" << endl;
					else acc << arrays.whole[i] << '\t' << whole  << "\tC" << endl;
#ifndef SINGLE
					pthread_mutex_lock(&tablelock);
					totals_info *info=&table[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
					pthread_mutex_unlock(&tablelock);
				}
				pthread_mutex_lock(&outlock);
				outfile << acc.str();
				pthread_mutex_unlock(&outlock);
#else
					totals_info *info=&table[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
				}
				outfile << acc.str();
#endif
			}
		}
	}
	return;
}

void query_sf(string &line) {																								//perform intersection for strand-aware flagged matching
	string desc1;
	double desc2;
	string chr;
	long physical_start;
	long physical_end;
	string strand;
	string whole;
	bool qstartcmp;
	bool qendcmp;
	ptr_entry arrays;
	ostringstream acc;
	istringstream temp1(line);
	temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end >> strand;
	if(temp1.fail()) {
		cout << "Query File contains bad line, skipping: " << line << endl;
	}
	else {
		ostringstream temp2;
		temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end << '\t' << strand;
		whole=temp2.str();
		if(db_s.find(chr)!=db_s.end()) {
			size_t max=db_s[chr].dbwhole.size();
			arrays.desc1=&db_s[chr].dbdesc1[0];
			arrays.desc2=&db_s[chr].dbdesc2[0];
			arrays.physical_start=&db_s[chr].dbphysical_start[0];
			arrays.physical_end=&db_s[chr].dbphysical_end[0];
			arrays.strand=&db_s[chr].dbstrand[0];
			arrays.whole=&db_s[chr].dbwhole[0];
			for(size_t i=0;i<max;i++) {
				if(arrays.physical_end[i]<physical_start || arrays.physical_start[i]>physical_end) continue;
				qstartcmp=(physical_start>=arrays.physical_start[i]);
				qendcmp=(physical_end<=arrays.physical_end[i]);
				if(qstartcmp && qendcmp) acc << arrays.whole[i] << '\t' << whole  << "\tB";
				else if(qstartcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tS";
				else if(qendcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tE";
				else acc << arrays.whole[i] << '\t' << whole  << "\tC";
				if(strand==arrays.strand[i]) {																				//flag as sense or anti-sense
					acc << "\tS\n";
				}
				else {
					acc << "\tA\n";
				}
#ifndef SINGLE
				pthread_mutex_lock(&tablelock);
				totals_info *info=&table[arrays.desc1[i]];
				info->total+=desc2;
				if(desc2>info->most) {
					info->most=desc2;
					info->most_loc.clear();
					info->most_loc.push_back(physical_start);
				}
				else if(desc2==info->most) {
					info->most_loc.push_back(physical_start);
				}
				pthread_mutex_unlock(&tablelock);
			}
			pthread_mutex_lock(&outlock);
			outfile << acc.str();
			pthread_mutex_unlock(&outlock);
#else
				totals_info *info=&table[arrays.desc1[i]];
				info->total+=desc2;
				if(desc2>info->most) {
					info->most=desc2;
					info->most_loc.clear();
					info->most_loc.push_back(physical_start);
				}
				else if(desc2==info->most) {
					info->most_loc.push_back(physical_start);
				}
			}
			outfile << acc.str();
#endif
		}
	}
	return;
}

void query_ss(string &line) {																								//perform intersection for strand-aware separate-output matching
	string desc1;
	double desc2;
	string chr;
	long physical_start;
	long physical_end;
	string strand;
	string whole;
	bool qstartcmp;
	bool qendcmp;
	ptr_entry arrays;
	ostringstream acc;
	ostringstream acc2;
	istringstream temp1(line);
	temp1 >> desc1 >> desc2 >> chr >> physical_start >> physical_end >> strand;
	if(temp1.fail()) {
		cout << "Query File contains bad line, skipping: " << line << endl;
	}
	else {
		ostringstream temp2;
		temp2 << desc1 << '\t' << desc2 << '\t' << chr << '\t' << physical_start << '\t' << physical_end << '\t' << strand;
		whole=temp2.str();
		if(db_s.find(chr)!=db_s.end()) {
			size_t max=db_s[chr].dbwhole.size();
			arrays.desc1=&db_s[chr].dbdesc1[0];
			arrays.desc2=&db_s[chr].dbdesc2[0];
			arrays.physical_start=&db_s[chr].dbphysical_start[0];
			arrays.physical_end=&db_s[chr].dbphysical_end[0];
			arrays.strand=&db_s[chr].dbstrand[0];
			arrays.whole=&db_s[chr].dbwhole[0];
			for(size_t i=0;i<max;i++) {
				if(arrays.physical_end[i]<physical_start || arrays.physical_start[i]>physical_end) continue;
				qstartcmp=(physical_start>=arrays.physical_start[i]);
				qendcmp=(physical_end<=arrays.physical_end[i]);
				if(strand==arrays.strand[i]) {
					if(qstartcmp && qendcmp) acc << arrays.whole[i] << '\t' << whole  << "\tB" << endl;
					else if(qstartcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tS" << endl;
					else if(qendcmp == 1) acc << arrays.whole[i] << '\t' << whole  << "\tE" << endl;
					else acc << arrays.whole[i] << '\t' << whole  << "\tC" << endl;
#ifndef SINGLE
					pthread_mutex_lock(&tablelock);
					totals_info *info=&table[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
					pthread_mutex_unlock(&tablelock);
				}
				else {
					if(qstartcmp && qendcmp) acc2 << arrays.whole[i] << '\t' << whole  << "\tB" << endl;
					else if(qstartcmp == 1) acc2 << arrays.whole[i] << '\t' << whole  << "\tS" << endl;
					else if(qendcmp == 1) acc2 << arrays.whole[i] << '\t' << whole  << "\tE" << endl;
					else acc2 << arrays.whole[i] << '\t' << whole  << "\tC" << endl;
					pthread_mutex_lock(&tablelock2);
					totals_info *info=&table2[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
					pthread_mutex_unlock(&tablelock2);
				}
			}
			pthread_mutex_lock(&outlock);
			outfile << acc.str();
			outfile2 << acc2.str();
			pthread_mutex_unlock(&outlock);
#else
					totals_info *info=&table[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
				}
				else {
					if(qstartcmp && qendcmp) acc2 << arrays.whole[i] << '\t' << whole  << "\tB" << endl;
					else if(qstartcmp == 1) acc2 << arrays.whole[i] << '\t' << whole  << "\tS" << endl;
					else if(qendcmp == 1) acc2 << arrays.whole[i] << '\t' << whole  << "\tE" << endl;
					else acc2 << arrays.whole[i] << '\t' << whole  << "\tC" << endl;
					totals_info *info=&table2[arrays.desc1[i]];
					info->total+=desc2;
					if(desc2>info->most) {
						info->most=desc2;
						info->most_loc.clear();
						info->most_loc.push_back(physical_start);
					}
					else if(desc2==info->most) {
						info->most_loc.push_back(physical_start);
					}
				}
			}
			outfile << acc.str();
			outfile2 << acc2.str();
#endif
		}
	}
	return;
}

#ifndef SINGLE
void *t_query(void *file) {												//called for every thread requested by user, read single lines from query file, pass to appropriate intersection function
	ifstream *qfile=reinterpret_cast<ifstream*>(file);
	string line;
	switch(s_val) {
	case 0:
		pthread_mutex_lock(&linelock);
		getline(*qfile,line);
		while(!qfile->eof()) {
			pthread_mutex_unlock(&linelock);
			query_i(line);
			pthread_mutex_lock(&linelock);
			getline(*qfile,line);
		}
		pthread_mutex_unlock(&linelock);
		break;
	case 1:
	case 2:
		pthread_mutex_lock(&linelock);
		getline(*qfile,line);
		while(!qfile->eof()) {
			pthread_mutex_unlock(&linelock);
			query_split(line);
			pthread_mutex_lock(&linelock);
			getline(*qfile,line);
		}
		pthread_mutex_unlock(&linelock);
		break;
	case 3:
		pthread_mutex_lock(&linelock);
		getline(*qfile,line);
		while(!qfile->eof()) {
			pthread_mutex_unlock(&linelock);
			query_sf(line);
			pthread_mutex_lock(&linelock);
			getline(*qfile,line);
		}
		pthread_mutex_unlock(&linelock);
		break;
	case 4:
		pthread_mutex_lock(&linelock);
		getline(*qfile,line);
		while(!qfile->eof()) {
			pthread_mutex_unlock(&linelock);
			query_ss(line);
			pthread_mutex_lock(&linelock);
			getline(*qfile,line);
		}
		pthread_mutex_unlock(&linelock);
		break;
	}
	pthread_exit(NULL);
}
#endif

int main(int argc,char** args) {
	ifstream infile,qfile;
	ofstream totalfile,totalfile2;
	string line;
	string dbf,qf,of,of2;
	int opt,dummy;
	int t=1;
	istringstream temp;
	while((opt=getopt_long(argc,args,"s:h:t:",long_options,&dummy))!=-1) {
		switch(opt) {
		case 'h':
			usage();
			return(0);
		case 's':
			if(!strcmp(optarg,"i")) {
				s_val=0;
			}
			else if(!strcmp(optarg,"s")) {
				s_val=1;
			}
			else if(!strcmp(optarg,"o")) {
				s_val=2;
			}
			else if(!strcmp(optarg,"bf")) {
				s_val=3;
			}
			else if(!strcmp(optarg,"bs")) {
				s_val=4;
			}
			else {
				cout << "Error: \"" << optarg << "\" is not a supported argument for the \"-s [ --strands ]\" option\n";
				usage();
				return(1);
			}
			break;
		case 't':
			temp.str(optarg);
			temp >> t;
			if(temp.fail()) {
				cout << "Error: -t argument must be an integer value\n";
				usage();
				return(1);
			}
			break;
		default:
			usage();
			return(1);
		}
	}
	dummy=argc-optind;
	if(dummy==0) {
		cout << "Error: DB file name must be specified\n";
		usage();
		return(1);
	}
	else if(dummy==1) {
		cout << "Error: query file name must be specified\n";
		usage();
		return(1);
	}
	else if(dummy==2) {
		cout << "Error: output file name must be specified\n";
		usage();
		return(1);
	}
	dbf=args[optind];
	qf=args[optind+1];
	of=args[optind+2];
	infile.open(dbf.c_str());
	if(infile.fail()) {
		cout << "Error: Could not open DB file \"" << dbf << "\"\n";
		return(1);
	}
	qfile.open(qf.c_str());
	if(qfile.fail()) {
		cout << "Error: Could not open query file \"" << qf << "\"\n";
		return(1);
	}
	switch(s_val) {																			//store contents of DB file using appropriate caching function
	case 0:
		cache_i(infile);
		break;
	case 1:
	case 2:
		cache_split(infile);
		break;
	default:
		cache_s(infile);
		break;
	}
	if(s_val!=0) {
		switch(strand_list.size()) {														//check strands listed in DB file
		case 0:
			cerr << "Error: DB File does not contain a strand identifier column" << endl;
			return(1);
		case 1:
			strand_list.insert("dummy");													//create dummy entry if DB file only contains single strand
			break;
		case 2:
			break;
		default:
			cerr << "Error: DB File contains more than two strand identifiers" << endl;
			return(1);
		}
		switch(s_val) {
		case 1:
			strand_map[*strand_list.begin()]=*strand_list.begin();							//for same-strand matching
			strand_map[*strand_list.rbegin()]=*strand_list.rbegin();
			break;
		case 2:
			strand_map[*strand_list.begin()]=*strand_list.rbegin();							//for opposite-strand matching
			strand_map[*strand_list.rbegin()]=*strand_list.begin();
			break;
		case 4:
			of2=of+"_antisense";															//create antisense output file for separate-output matching
			outfile2.open(of2.c_str());
			if(outfile2.fail()) {
				cout << "Error: Could not create output file \"" << of << "\"\n";
				return(1);
			}
			of+="_sense";
			break;
		}
	}
	outfile.open(of.c_str());
	if(outfile.fail()) {
		cout << "Error: Could not create output file \"" << of << "\"\n";
		return(1);
	}
#ifndef SINGLE
	if(t!=1) {
		pthread_t tid[t];
		int ti;
		for(ti=0;ti<t;ti++) {															//create number of threads specified by user, each reads a single line at a time, performs the appropriate intersection
			pthread_create(&tid[ti],NULL,t_query,reinterpret_cast<void*>(&qfile));
		}
		for(ti=0;ti<t;ti++) {
			pthread_join(tid[ti],NULL);													//all threads exit when no lines remain to be read from the query file
		}
	}
	else {
		getline(qfile,line);
		while(!qfile.eof()) {															//for single-thread runs, read a line at a time from the query file, pass to appropriate intersection function
			switch(s_val) {
			case 0:
				query_i(line);
				break;
			case 1:
			case 2:
				query_split(line);
				break;
			case 3:
				query_sf(line);
				break;
			case 4:
				query_ss(line);
				break;
			}
			getline(qfile,line);
		}
	}
#else
	getline(qfile,line);
	while(!qfile.eof()) {
		switch(s_val) {
		case 0:
			query_i(line);
			break;
		case 1:
		case 2:
			query_split(line);
			break;
		case 3:
			query_sf(line);
			break;
		case 4:
			query_ss(line);
			break;
		}
	getline(qfile,line);
	}
#endif
	string f;
	f=of+"_total";
	totalfile.open(f.c_str());
	if(totalfile.fail()) {
		cout << "Error: Could not create output file \"" << f << "\"\n";
		return(1);
	}
	for(map<string,totals_info>::iterator i=table.begin();i!=table.end();i++) {												//print results to *_total files, separate tied most hit locations with colons
		totalfile << i->first << '\t' << i->second.total << '\t' << i->second.most << '\t' << i->second.chr << '\t';
		for(vector<long>::iterator j=i->second.most_loc.begin();j!=i->second.most_loc.end();j++) {
			if(j!=i->second.most_loc.begin()) {
				totalfile << ":" << *j;
			}
			else {
				totalfile << *j;
			}
		}
		totalfile << endl;
	}
	if(s_val==4) {
		f=of2+"_total";
		totalfile2.open(f.c_str());
		if(totalfile2.fail()) {
			cout << "Error: Could not create output file \"" << f << "\"\n";
			return(1);
		}
		for(map<string,totals_info>::iterator i=table2.begin();i!=table2.end();i++) {
				totalfile2 << i->first << '\t' << i->second.total << '\t' << i->second.most << '\t' << i->second.chr << '\t';
				for(vector<long>::iterator j=i->second.most_loc.begin();j!=i->second.most_loc.end();j++) {
					if(j!=i->second.most_loc.begin()) {
						totalfile2 << ":" << *j;
					}
					else {
						totalfile2 << *j;
					}
				}
				totalfile2 << endl;
		}
		totalfile2.close();
	}
	totalfile.close();
	qfile.close();
	outfile.close();
	return(0);
}
