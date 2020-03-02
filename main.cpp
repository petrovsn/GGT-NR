#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <map>
#include <set>
#include <math.h>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include <sstream>
#include <thread>
#include <future>


#include <ctime>
#include <cstdio>

#include<Graph2.h>
#include <Interface.h>

using namespace std;
int main(int argc, const char* argv[])
{
	Graph GMAIN;
	//GMAIN.LoadFromGFA("data/001_best_spades_graph.gfa");
	GMAIN.LoadFromGFA("data/assembly_graph_with_scaffolds.gfa");
	//GMAIN.LoadFromGFA("data/asm_nature.gfa");



	Node tmp = GMAIN.Body[-2428];
	GMAIN.BuildIndex(5, 13);
	vector<pair<string, string>> reads = loadmultifasta2("data/testEcSim.fa");
	string read = "CGGTGGCGGGAATAACGGTGACCTTCACCATGCCACAGGACGTGGCGGCAAACTTTACCCTCGAAAATAACGGTATTGCCATCACCCAGGCCAATGGGGAAGCGCATGTCACGCTCAAAGGTAAAAAAGCGGGCACGCATACGGTTACCGTTTTTTTTTTTTTTTT";
	string revread = reverse(read);
	FAlignment tmp3 = GMAIN.AlignHashSmWtmn(read);
	FAlignment tmp4 = GMAIN.AlignHashSmWtmn(revread);
	cout << "loaded";
}

/*int main(int argc, const char* argv[])
{
	string REF_IN = argv[1];
	string READS_IN = argv[2];
	int THREADS = stoi(argv[3]);


	//string REF_IN = "data/Ecoli_O157.fasta";
	//string READS_IN = "data/sample.fa";
	//int THREADS = 20;

	Graph GMAIN;
	//GMAIN.LoadFromGFA("data/assembly.gfa");

	GMAIN.LoadFromGFA(REF_IN);
	GMAIN.BuildIndex(5, 13);

	Results res = Align2GMAIN(GMAIN, READS_IN, 0.85, THREADS);

	cout << "Aligned:" << res.sumAlignment.size();

	for (auto aln : res.sumAlignment)
	{
		GMAIN.GetCoverage(aln.second);
	}

	ofstream var = ofstream("var2.txt");

	for (auto t : res.VarMap)
	{
		//GVariation g = filtered_vars[i];
		GVariation g = t.first;
		var << g.ID1 << '.' << g.pos1 << '.' << g.ID2 << '.' << g.pos2 << '\t' << g.alt << '\t';
		var << res.VarMapStringCount[g].size() << '\n';
		for (int i = 0; i < res.VarMapStringCount[g].size(); i++)
		{
			var << '\t' << res.VarMapStringCount[g][i] << endl;
		}
	}
	var.close();

	vector<GVariation> filtered_vars = GMAIN.FilterVarMap(res.VarMap);

	var = ofstream("var2_filtered.txt");

	for (auto t : res.VarMap)
	{
		//GVariation g = filtered_vars[i];
		GVariation g = t.first;
		var << g.ID1 << '.' << g.pos1 << '.' << g.ID2 << '.' << g.pos2 << '\t' << g.alt << '\t';
		var << res.VarMapStringCount[g].size() << '\n';
		for (int i = 0; i < res.VarMapStringCount[g].size(); i++)
		{
			var << '\t' << res.VarMapStringCount[g][i] << endl;
		}
	}
	var.close();


	//GMAIN.InitOld2NewIds();
	//GMAIN.UpgradeGraph(filtered_vars);
	//GMAIN.Save("newgraph.mygraph.txt");
	//GMAIN.ReInitGraph();
	//GMAIN.BuildIndex(5, 13);
}

int main(int argc, const char* argv[])
{


	//vector<pair<string, string>> reads = loadmultifasta2("data/outSRR65fs.bam.fa");
	vector<pair<string, string>> reads = loadmultifasta2("data/testEcSim.fa");
	//vector<pair<string, string>> reads = loadmultifasta2("data/SRR65filtered_sorted_2879898.fa");




	Graph GMAIN = Graph();
	//GMAIN.LoadFromGFA("data/001_best_spades_graph.gfa");
//	GMAIN.LoadFromGFA("data/assembly_graph_with_scaffolds.gfa");
	GMAIN.LoadFromGFA("data/asm_nature.gfa");
	GMAIN.BuildIndex(5, 13);

	map<string, vector<FAlignment>> results;

	clock_t res = 0;
	float treshold = 0.8;

	int aligned = 0;
	int nonaligned = 0;

	clock_t time_avg;

	//ofstream fout = ofstream("data/nonaligned5.fa");

	unordered_map<GVariation, int, MyHashFunction> VarMap;
	map<string, FAlignment> sumAlignment;
	int maxlen = reads.size();
	for (int i = 39014; i < reads.size(); i++)
	{
		cout << "read#: " << i << '\r';
		clock_t time_a = clock();
		string forread = reads[i].second;
		string revread = reverse(reads[i].second);
		FAlignment aln = GMAIN.AlignHashSmWtmn(forread);
		FAlignment aln_rev = GMAIN.AlignHashSmWtmn(revread);

		if (aln_rev.score > aln.score)
		{
			aln = aln_rev;
		}

		if (aln.tresh > treshold)
		{
			for (int j = 0; j < aln.vars.size(); j++)
			{
				GVariation key = aln.vars[j];
				VarMap[key] = VarMap[key] + 1;
			}

			sumAlignment[reads[i].first] = aln;
			aligned = aligned + 1;
		}
		else
		{
			nonaligned = nonaligned + 1;
			//fout << '>' << reads[i].first << '\n';
			//fout << reads[i].second << '\n';
		}


		clock_t time_b = clock() - time_a;
		res = res + time_b;
		//cout << res * 1.0 / (i + 1) << '\t' << i <<'/'<< maxlen << '\r';
		time_avg = res * 1.0 / (i + 1);
	}
	//	fout.close();
	cout << "Avgtime:\t" << time_avg << endl;
	cout << "Total:\t" << reads.size() << endl;
	cout << "aligned:\t" << aligned << endl;
	cout << "no aligned:\t" << nonaligned << endl;


	ofstream var = ofstream("data/var_local.txt");

	for (auto t : VarMap)
	{
		//GVariation g = filtered_vars[i];
		GVariation g = t.first;
		//var << g.ID1 << '.' << g.pos1 << '.' << g.ID2 << '.' << g.pos2 << '\t' << g.alt << '\n';
	}
	var.close();




	//GMAIN.ReInitGraph();
	//GMAIN.BuildIndex(5, 13);
}

/*int main(int argc, const char* argv[])
{

	Graph GMAIN;

	GMAIN.BuildIndex(5, 13);
	GMAIN.LoadFromGFA("data/assembly_graph_with_scaffolds.gfa");


	Results res = Align2GMAIN(GMAIN, "data/testEcSim.fa", 0.85, 4);

	cout << "Aligned:" << res.sumAlignment.size();

/*	for (auto aln : res.sumAlignment)
	{
		GMAIN.GetCoverage(aln.second);
	}


}

int main(int argc, const char* argv[])
{
	Graph GMAIN;
	string REF_IN = "data/Ecoli_O157.fasta";
	string READS_IN = "data/SRR65filtered_sorted_pos3434_fixed.fa";
	int THREADS = 6;


	string ref = loadfasta(REF_IN);

	GMAIN.LoadReference(ref);

	int step = 1000;
	vector<GVariation> gvars;
	for (int i = 2; i < ref.size(); i += 20)
	{
		GVariation tmpgvar;
		tmpgvar.ID1 = 0;
		tmpgvar.ID2 = 0;
		tmpgvar.pos1 = i;
		tmpgvar.pos2 = i + 2;
		tmpgvar.alt = "a";

		gvars.push_back(tmpgvar);
	}

	GMAIN.InitOld2NewIds();
	GMAIN.UpgradeGraph(gvars);

	vector<pair<string, string>> reads = loadmultifasta2("data/SRR65filtered_sorted_pos3434_fixed.fa");

	GMAIN.BuildIndex(5, 13);

 	float treshold = 0.8;

	int aligned = 0;
	int nonaligned = 0;

	clock_t time_avg;

	//ofstream fout = ofstream("data/nonaligned5.fa");
	clock_t res = 0;

	unordered_map<GVariation, int, MyHashFunction> VarMap;
	map<string, FAlignment> sumAlignment;
	int maxlen = reads.size();
	for (int i = 0; i < reads.size(); i++)
	{
		clock_t time_a = clock();
		string forread = reads[i].second;
		string revread = reverse(reads[i].second);
		FAlignment aln = GMAIN.AlignHashSmWtmn(forread);
		FAlignment aln_rev = GMAIN.AlignHashSmWtmn(revread);

		clock_t time_b = clock() - time_a;
		res = res + time_b;
		time_avg = res * 1.0 / (i + 1);
	}

	cout << "Avgtime:\t" << time_avg*1.0/2 << endl;


	cout << "endl";
}

/*int main(int argc, const char* argv[])
{
	string REF_IN = ""; 
	string READS_IN = "";  
	int THREADS = 0; 

	bool verbose = false;
	float treshhold = 0.8;

	for (int i = 1; i < argc; i++)
	{
		cout << i << '\t' << argv[i] << endl;
		string str_argi = argv[i];
		if (str_argi == "-g")
		{
			REF_IN = argv[i + 1];
		}

		if (str_argi == "-f")
		{
			READS_IN = argv[i + 1];
		}
		
		if (str_argi == "-t")
		{
			THREADS = stoi(argv[i+1]);
		}

		if (str_argi == "-v")
		{
			verbose = true;
		}

		if (str_argi == "-c")
		{
			treshhold = stof(argv[i + 1]);
		}
	}

	cout << "Ref:   " << REF_IN << endl;
	cout << "Reads: " << READS_IN << endl;
	cout << "Theads:" << THREADS << endl;
	cout << "Tresh: " << treshhold << endl;
	strcmp("AA", "VVV");

	Graph GMAIN;

	GMAIN.LoadFromGFA(REF_IN);
	GMAIN.BuildIndex(5, 13);

	Results res = Align2GMAIN(GMAIN, READS_IN, treshhold, THREADS, verbose);

	cout << "Aligned: " << res.sumAlignment.size() << endl;
	cout << "Avgtime: " << res.time_avg << endl;
}
*/