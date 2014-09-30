// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <algorithm>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "evaluation.h"
#include <iostream>
#include <lemon/list_graph.h>
//#include <lemon/lp_base.h>
#include <lemon/lp.h>
#include <lemon/maps.h>
#include <mflib/MethylRead.hpp>
//#include <lemon/bipartite_matching.h>
//#include <lemon/concepts/bpgraph.h>
//lemon::Lp lp;


typedef lemon::ListDigraph Graph;
typedef Graph::Node RedNode;
typedef Graph::Node BlueNode;
typedef Graph::NodeIt RedNodeIt;
typedef Graph::NodeIt BlueNodeIt;
typedef Graph::NodeIt NodeIt;
typedef Graph::Node Node;
typedef Graph::Arc Edge;
typedef Graph::ArcIt EdgeIt;
typedef Graph::InArcIt InEdgeIt;
typedef Graph::ArcMap<float> LengthMap;

using namespace methylFlow;

std::map<RedNode, MethylRead*> read_map;
std::map<RedNode, float> abundance_map;

std::ofstream evalFile;
std::ofstream MCFFile;



using lemon::INVALID;

Graph g;
LengthMap length(g);

int chr;
float var;
float minCostFlowErr;

std::ifstream truePatternFile, estimatedPatternFile, mFile, shortReadFile;
std::ofstream weightFile, matchFile, matchAppFile, mEstimatedPercentageFile, mTruePercentageFile, mReadPercentageFile;
std::vector<MethylRead*> trueMethylData,estimatedMethylData, readMethylData;
std::vector<float> trueAbundanceData, estimatedAbundanceData;


/////Variables for reading MatchMatrix
int truePatternNum, estimatedPatternNum;
vector<int> abdncTrue, abdncEstimated, idTrue, idEstimated;
std::map<int, float> abdnc_map;
std::map<int, int> matchTrue_map;
std::map<int, int> matchEstimated_map;
std::map<int, float> weight_map;


//Varibales for computing the avarage error matrices

float methylCallError = 0;
float abndncError = 0;
int TP, FP, FN = 0;

vector<float> thr;
std::map<float, vector<float> > abdncErr_avg_map;
std::map<float, vector<float> > methylErr_avg_map;
std::map<float, vector<int> > TP_avg_map;
std::map<float, vector<int> > FN_avg_map;
std::map<float, vector<int> > FP_avg_map;


vector<int> cpgPos;
std::map<int, float> true_U_map, true_M_map, estimated_U_map, estimated_M_map, read_U_map, read_M_map;


Graph::Node addNode_Read(MethylRead *read, float abndnc)
{
    Graph::Node n = g.addNode();
    read_map[n] = read;
    abundance_map[n] = abndnc;
    return n;
}
void readTruePattern(int start, int end){
	// "s" is for start position, "e" is for end position
	
	string chr, methylString;
    string dummyLine;
    getline(truePatternFile, dummyLine);
	int i = 0;
	while(!truePatternFile.eof()) {
        int  s = -100, e = -100, cid, pid;
        float abundance = 0 ;

		i++;
		truePatternFile >> chr >> s >> e >> cid >> pid >> abundance >> methylString;
        //cout << "Ture read : "<< chr << " " << s << " " << e <<endl;
		if( s >= start && e <= end){
            MethylRead* m = new MethylRead(s, e-s+1);
            m->parseMethyl(methylString);
            
            //m->write();
            
            trueMethylData.push_back(m);
            trueAbundanceData.push_back(abundance);
		}
	}
    float sum =0;
    for(unsigned int j= 0; j<trueAbundanceData.size(); j++){
        sum += trueAbundanceData.at(j);
    }
    for(unsigned int j= 0; j<trueAbundanceData.size(); j++){
        trueAbundanceData.at(j) = trueAbundanceData.at(j)/sum ;
    }
    
}

void readEstimatedPattern(int start, int end){
	// "s" is for start position, "e" is for end position
	string  chr, methylString;
    string dummyLine;
    getline(estimatedPatternFile, dummyLine);
	int i = 0;
	while(!estimatedPatternFile.eof()) {
        int  s = -100, e = -100, cid, pid;
        float abundance = 0 ;

		i++;
		estimatedPatternFile >> chr >> s >> e >> cid >> pid >> abundance >> methylString;
        
        //cout << "Estimated read : "<< chr << " " << s << " " << e << "  " << methylString << endl;
		if( s >= start -1  &&  e <= end){
            
            MethylRead* m = new MethylRead(s, e-s+1);
            //cerr << "befor parse" << endl;
            
            m->parseMethyl(methylString);
            //cerr << "befor write" << endl;
            //m->write();
            //cerr << "after write" << endl;
            
            estimatedMethylData.push_back(m);
            estimatedAbundanceData.push_back(abundance);
		}
	}
    float sum =0;
    for(unsigned int j= 0; j<estimatedAbundanceData.size(); j++){
        sum += estimatedAbundanceData.at(j);
    }
    for(unsigned int j= 0; j<estimatedAbundanceData.size(); j++){
        estimatedAbundanceData.at(j) = estimatedAbundanceData.at(j)/sum ;
    }
    
}

void readShortRead(int start, int end){
    int  s, l;
	string readId, methylString, strand, etc;
    //string dummyLine;
    //getline(truePatternFile, dummyLine);
	int i = 0;
	while(!shortReadFile.eof()) {
		i++;
		shortReadFile >> readId >> s >> l >> strand >> methylString >> etc;
        //cout << "Ture read : "<< chr << " " << s << " " << e <<endl;
		if( s >= start - 1 && s + l -1  <= end){
            MethylRead* m = new MethylRead(s, s + l -1);
            m->parseMethyl(methylString);
            
            //m->write();
            
            readMethylData.push_back(m);
		}
	}
}

float cost(Graph::Node u, Graph::Node v) {
    
    
    MethylRead* readU = read_map[u];
    MethylRead* readV = read_map[v];
    
    
    int common = 0 ;
    //cout << "Distance start" << endl;
    int match =  readU->distance(readV, common);
    //readU->write();
    //readV->write();
    //cout <<  "node " << g.id(u) << ", " << g.id(v) << endl;
    //cout << "match " << match <<  endl;
    
    //cout << "common " << common <<  endl;
    
    int mismatch = readU->cpgOffset.size() + readV->cpgOffset.size() - match - common ;
    //int mismatch = common - match;
    int totalCpG = readU->cpgOffset.size() + readV->cpgOffset.size() - common;
    //cout << "mismatch " << mismatch <<  endl;
    
    //cout << "totalCpG " << totalCpG <<  endl;
    
    
    //cout << "cost " << (float(mismatch) / totalCpG) <<  endl;
    
    //return (float(mismatch) / common) ;
    
    return (float(mismatch) / totalCpG) ;
}

void buildGraph() {
    
    //cout << "trueMethylSize = " << trueMethylData.size() << endl;
    //cout << "estimatedMethylSize = " << estimatedMethylData.size() << endl;
    for (int i=0; i < trueMethylData.size(); i++) {
       // cout << "Add true node " << i << endl;
        addNode_Read(trueMethylData.at(i), trueAbundanceData.at(i));
    }
    
    for (int i=0; i < estimatedMethylData.size(); i++) {
       // cout << "Add estimated node " << i << endl;
        addNode_Read(estimatedMethylData.at(i), estimatedAbundanceData.at(i));
    }
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        for (RedNodeIt v(g); v != INVALID; ++v) {
            if (g.id(u) < trueMethylData.size() && g.id(v) >= trueMethylData.size()) {
                //cout << "Add arc " << g.id(u) << " " <<g.id(v) << endl;
                
                g.addArc(u, v);
            }
        }
    }
    
    cout << "set weight" << endl;
    
    ///////////Set Weights///////////////
    
    for (EdgeIt e(g); e != INVALID; ++e) {
        RedNode u = g.source(e);
        BlueNode v = g.target(e);
        length[e] = cost(u, v);
        
        //cout << "Cost arc (" << g.id(u) <<", " << g.id(v) << ") = " << length[e] << endl;
    }
    
    cout << "build done" << endl;
}

void writeMatchMatrix() {
    // weightFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/weight.txt");
    // matchFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/match.txt");
    
    weightFile << trueMethylData.size() << "\t" << estimatedMethylData.size() << endl;
    matchFile << trueMethylData.size() << "\t" << estimatedMethylData.size() << endl;
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        if (g.id(u) < trueMethylData.size()) {
            weightFile << g.id(u) << "\t" << abundance_map[u] << endl;;
            matchFile << g.id(u) << "\t" << abundance_map[u] << endl;
        }
    }
    
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        if (g.id(u) >= trueMethylData.size()) {
            weightFile << g.id(u) << "\t" << abundance_map[u] << endl;;
            matchFile << g.id(u) << "\t" << abundance_map[u] << endl;
        }
    }
    
    
    
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        
        if (g.id(u) < trueMethylData.size()) {
            double min = 1.001;
            int argmin = -1;
            
            
            for (RedNodeIt v(g); v != INVALID; ++v) {
                if (g.id(v) >= trueMethylData.size()) {
                    float weight = cost(u, v);
                    weightFile << var << "\t" << g.id(u) << "\t" << g.id(v) << "\t" << weight << endl;
                    
                    if(min > weight ){
                        min = cost(u,v);
                        argmin = g.id(v);
                    }
                }
            }
            
            matchFile << var << "\t" << g.id(u) << "\t" << argmin  << "\t" << min << endl;
            matchAppFile << var << "\t" << g.id(u) << "\t" << argmin  << "\t" << min << endl;

        }
    }
    matchFile.close();
    weightFile.close();
    matchAppFile.close();
    
}


void readMatchMatrix() {
    //reading the match file for pruning the matched edges
    int id, matchId;
    float abdnc, weight, var;
    
    // mFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/match.txt");
    mFile >> truePatternNum  >> estimatedPatternNum ;
    for (int i = 0; i < truePatternNum ; i++){
        mFile >> id >> abdnc;
        idTrue.push_back(id);
        abdnc_map[id] = abdnc;
    }
    for (int i = 0; i < estimatedPatternNum ; i++){
        mFile >> id >> abdnc;
        idEstimated.push_back(id);
        abdnc_map[id] = abdnc;
    }
    
    for (int i=0; i< truePatternNum ; i++ ){
        mFile >> var >> id >> matchId >> weight ;
        matchTrue_map[id] = matchId;
        weight_map[id] = weight;
        matchEstimated_map[matchId] = id;
    }
    
    mFile.close();
}

void computeErrorMatrix(double threshold) {
    float methylCallError = 0;
    float abndncError = 0;
    
    int match = 0;
    std::set<int> matchSet;
    for (int i = 0 ; i < truePatternNum ; i ++){
        if(weight_map[i] < threshold) {
            abndncError += pow(double(abdnc_map[i] - abdnc_map[matchTrue_map[i]]),2)/1.0;
            matchSet.insert(matchTrue_map[i]);
            methylCallError += weight_map[i];
            match++;
        }
        /*else{
         abndncError += pow(double(abdnc_map[i]),2) /10000;
         abndncError += pow(double(abdnc_map[matchTrue_map[i]]),2) /10000;
         //methylCallError += weight_map[i];
         }*/
    }
    
    for(int i = truePatternNum; i < truePatternNum + estimatedPatternNum ; i++){
        if (matchEstimated_map.find(i) == matchEstimated_map.end()) {
            abndncError += pow(double(abdnc_map[i]),2) /1.0;
        }
        
    }
    // for case with less estimated pattern than true pattern
    
    
    int TP = match;
    int FN = truePatternNum - match;
    int FP = estimatedPatternNum - matchSet.size();
    
    //cout << "match = " << match << endl;
    //cout << "FP = " << FP << endl;
    //cout << "estimatedPatternNum = " << estimatedPatternNum << endl;
    //cout << "truePatternNum = " << truePatternNum << endl;

    methylCallError = methylCallError/match;
    abndncError = sqrt(abndncError/TP);
    if( !TP == 0 ){
        evalFile << threshold << "\t" << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
        thr.push_back(threshold);
        abdncErr_avg_map[threshold].push_back(abndncError);
        methylErr_avg_map[threshold].push_back(methylCallError);
        TP_avg_map[threshold].push_back(TP);
        FN_avg_map[threshold].push_back(FN);
        FP_avg_map[threshold].push_back(FP);
    }
    else{
        evalFile << threshold << "\t" << 0 << "\t" << 0 << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
        thr.push_back(threshold);
        abdncErr_avg_map[threshold].push_back(0);
        methylErr_avg_map[threshold].push_back(0);
        TP_avg_map[threshold].push_back(TP);
        FN_avg_map[threshold].push_back(FN);
        FP_avg_map[threshold].push_back(FP);
    }
    
    //evalFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/eval.txt", std::ios_base::app);
    
    //evalFile << "abndncError " << "\t"  << "methylCall Error " << "\t" << "TP "<< "\t" << "FN " << "\t" << "FP " << endl;
    
    
    //std:cout << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
    
    
    
    
}

float computeMinCostFlowError(){
    if (estimatedMethylData.size() == 0) {
        return 1;
    }
    
    lemon::Lp lp;
    lp.min();
    
    std::map<Edge, lemon::Lp::Col> f;
    
    std::map<RedNode, lemon::Lp::Col> alpha;
    
    //  std::map<Graph::BlueNode, lemon::Lp::Col> beta;
    
    
    lemon::Lp::Expr obj;
    
    
    for (EdgeIt e(g); e != INVALID; ++e) {
        
        f[e] = lp.addCol();
        obj += f[e]* length[e];
        lp.colLowerBound(f[e], 0.0);
        //cerr << "length " <<  g.id(g.source(e)) <<"\t" << g.id(g.target(e)) << "\t" << length[e] << endl;
        
        
    }
    
    for (NodeIt u(g); u != INVALID; ++u) {
        //cerr <<  g.id(u) << "\t" << abundance_map[u] << endl;
        lemon::Lp::Expr c;
        if (g.id(u) < trueMethylData.size()) {
            for (lemon::ListDigraph::OutArcIt arcIt(g, u); arcIt != INVALID; ++arcIt) {
                
                Edge arc(arcIt);
                c += f[arc];
                
            }
            lp.addRow(c >= abundance_map[u]);
            
            
        }
        else{
            for (lemon::ListDigraph::InArcIt arcIt(g, u); arcIt != INVALID; ++arcIt) {
                
                Edge arc(arcIt);
                c += f[arc];
                
            }
            lp.addRow(c >= abundance_map[u]);
            
        }
    }
    
    
    lp.obj(obj);
    lp.solve();
    
    float v = lp.primal();
    //cout << "mincostflow error = " << v << endl;
    for (EdgeIt e(g); e != INVALID; ++e) {
        lemon::Lp::Col col = f[e];
        float val = lp.primal(col);
        //cerr << "primal " << g.id(g.source(e)) << "\t " << g.id(g.target(e)) << "\t" << val << endl;
    }
    
    return v;
    
    
    
}

void computeMethylPercentage(){
    for (unsigned int i=0; i < trueMethylData.size(); i++) {
        for (unsigned int j=0; j < trueMethylData.at(i)->cpgOffset.size(); j++) {
            int position = trueMethylData.at(i)->cpgOffset.at(j) + trueMethylData.at(i)->start();
            bool meth = trueMethylData.at(i)->methyl.at(j);
            if (find (cpgPos.begin(), cpgPos.end(), position) == cpgPos.end()) {
                cpgPos.push_back(position);
            }
            if (meth) {
                true_M_map[position] += trueAbundanceData.at(i);
                //cerr << trueAbundanceData.at(i) << "\t" << "true_M_map[position]= " << true_M_map[position] << endl;

            }
            else{
                true_U_map[position] += trueAbundanceData.at(i);
               // cerr << trueAbundanceData.at(i) << "\t" << "true_U_map[position]= " << true_U_map[position] << endl;

            }
        }
    }
    for (unsigned int i=0; i < estimatedMethylData.size(); i++) {
        //cerr << "estimated " << endl;
        for (unsigned int j=0; j < estimatedMethylData.at(i)->cpgOffset.size(); j++) {
            int position = estimatedMethylData.at(i)->cpgOffset.at(j) + estimatedMethylData.at(i)->start();
            bool meth = estimatedMethylData.at(i)->methyl.at(j);
            if (find (cpgPos.begin(), cpgPos.end(), position) == cpgPos.end()) {
                cpgPos.push_back(position);
            }
            if (meth) {
                estimated_M_map[position] += estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" << "estimated_M_map[position]= " << estimated_M_map[position] << endl;
            }
            else{
                estimated_U_map[position] += estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" <<  "estimated_U_map[position]= " << estimated_U_map[position] << endl;

            }
        }
    }
    
    for (unsigned int i=0; i < readMethylData.size(); i++) {
        //cerr << "short read " << endl;
        for (unsigned int j=0; j < readMethylData.at(i)->cpgOffset.size(); j++) {
            int position = readMethylData.at(i)->cpgOffset.at(j) + readMethylData.at(i)->start();
            bool meth = readMethylData.at(i)->methyl.at(j);
            if (find (cpgPos.begin(), cpgPos.end(), position) == cpgPos.end()) {
                cpgPos.push_back(position);
            }
            if (meth) {
                read_M_map[position] += 1;
                //cerr << "1" << "\t" << "read_M_map[position]= " << read_M_map[position] << endl;
            }
            else{
                read_U_map[position] += 1;
                //cerr << "1" << "\t" <<  "read_U_map[position]= " << read_U_map[position] << endl;
                
            }
        }
    }

    
}

void writeMethylPercentageMatrix(){
    mTruePercentageFile << "pos" << "\t" << "unmethyl" << "\t" << " methyl" << endl;
    for (unsigned int i = 0; i < cpgPos.size(); i++) {
        mTruePercentageFile << cpgPos.at(i) << "\t"  << true_U_map[cpgPos.at(i)] << "\t"  << true_M_map[cpgPos.at(i)] << endl;
    }
    
    
    
    mEstimatedPercentageFile << "pos" << "\t" << "unmethyl" << "\t" << " methyl" << endl;

    for (unsigned int i = 0; i < cpgPos.size(); i++) {
        mEstimatedPercentageFile << cpgPos.at(i) << "\t"  << estimated_U_map[cpgPos.at(i)] << "\t"  << estimated_M_map[cpgPos.at(i)] << endl;
    }
    
    
    mReadPercentageFile << "pos" << "\t" << "unmethyl" << "\t" << " methyl" << endl;

    for (unsigned int i = 0; i < cpgPos.size(); i++) {
        mReadPercentageFile << cpgPos.at(i) << "\t"  << read_U_map[cpgPos.at(i)] << "\t"  << read_M_map[cpgPos.at(i)] << endl;
    }
    
}




int main (int argc, char* argv[]) {
    
    
	int start, end;
    std::stringstream buffer;
    
	//cout << "start main" << endl;
    
	if(argc < 6){
		cout << "Please enter your input" << endl;
		return -1;
	}
	if (argc >= 6){
		start = atoi(argv[3]);
		end = atoi(argv[4]);
		var = atof(argv[5]);
        
	}
    
    std::string indirname = argv[1];
    
    buffer.str("");
    buffer << indirname << "/simPattern.txt";
    truePatternFile.open( buffer.str().c_str() );
    
    
    buffer.str("");
    buffer << indirname << "/patterns.tsv";
    estimatedPatternFile.open( buffer.str().c_str() );
    
    
    buffer.str("");
    buffer << indirname << "/shortRead.txt";
    shortReadFile.open( buffer.str().c_str() );
    
    
    std::string outdirname = argv[2];
    
    buffer.str("");
    buffer << outdirname << "/eval.txt";
    evalFile.open( buffer.str().c_str() ,std::ios_base::app);
    
    
    buffer.str("");
    buffer << outdirname << "/mcf.txt";
    MCFFile.open( buffer.str().c_str(), std::ios_base::app);
    
    
    buffer.str("");
    buffer << outdirname << "/match.txt";
    matchFile.open( buffer.str().c_str());
    
    buffer.str("");
    buffer << outdirname << "/matchApp.txt";
    matchAppFile.open( buffer.str().c_str(), std::ios_base::app);
    
    
    buffer.str("");
    buffer << outdirname << "/weight.txt";
    weightFile.open( buffer.str().c_str());
    
    
    buffer.str("");
    buffer << indirname << "/match.txt";
    mFile.open( buffer.str().c_str() );
    
    
    
    buffer.str("");
    buffer << outdirname << "/methylPercentageTrue.txt";
    mTruePercentageFile.open( buffer.str().c_str() );
    
    buffer.str("");
    buffer << outdirname << "/methylPercentageEstimated.txt";
    mEstimatedPercentageFile.open( buffer.str().c_str() );
    
    
    buffer.str("");
    buffer << outdirname << "/methylPercentageRead.txt";
    mReadPercentageFile.open( buffer.str().c_str() );
    
    
	//truePatternFile.open(argv[1]);
    //estimatedPatternFile.open(argv[2]);
    //evalFile.open(argv[3],std::ios_base::app);
    //MCFFile.open(argv[4],std::ios_base::app);
    
	
	//################     read the input .tsv data to the "line" number
	cout << "reading patterns " << start << endl;
	readTruePattern(start , end);
    readEstimatedPattern(start , end);
    cout << "reading short read " <<  endl;
    readShortRead(start, end);
    cout << "compute methyl Percentage  " << endl;

    computeMethylPercentage();
    cout << "write methyl Percentage  " << endl;
    writeMethylPercentageMatrix();
    
	cout << "build Graph " << start << endl;
    buildGraph();
    
    
    //########### write the weight of edges and corresponding matches for every true pattern ########
    /// the ffirst line is #truePattern , #estimated Pattern
    /// then each line is true Pattern id , abundance of pattern
    /// then each line is estimated pattern id , abundance of pattern
    // then the weight information of edges and matching information is seen in the rest of file
    
    minCostFlowErr =  computeMinCostFlowError();
    MCFFile << var << "\t" << minCostFlowErr << endl;
    cout << "min cost flow error " << minCostFlowErr << endl;
    
    
    writeMatchMatrix();
    
    
    readMatchMatrix();
    
    
    //// computing Error Metric ////////
    for (double thresh=0.1; thresh<0.4; thresh +=0.05)
        computeErrorMatrix(thresh);
    
    
    
    
    //}
    
    /*
     cout << "Solve LP " << start << endl;
     
     ////////////LP//////////////////////
     lemon::Lp lp;
     lp.max();
     
     //Graph::ArcMap<lemon::Lp::Col> x;
     std::map<Edge, lemon::Lp::Row> x;
     
     std::map<RedNode, lemon::Lp::Col> alpha;
     //  std::map<Graph::BlueNode, lemon::Lp::Col> beta;
     
     
     lemon::Lp::Expr obj;
     
     for (RedNodeIt r(g); r != INVALID; ++r) {
     alpha[r] = lp.addCol();
     //lp.colLowerBound(alpha[r], 0.0);
     obj += alpha[r];
     }
     
     
     
     for (EdgeIt e(g); e != INVALID; ++e) {
     
     //   x[e] = lp.addRow();
     RedNode r = g.source(e);
     BlueNode b = g.target(e);
     
     x[e] = lp.addRow(alpha[r] + alpha[b] - length[e] <= 0);
     }
     
     
     
     lp.obj(obj);
     lp.solve();
     
     float v = lp.primal();
     //cout << "dual : " << v << endl;
     
     float abndncError = 0;
     float methylCallError = 0;
     
     
     int match = 0;
     for (EdgeIt e(g); e != INVALID; ++e) {
     
     lemon::Lp::Row row = x[e];
     float dual = lp.dual(row);
     
     cout << " X(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << dual<< endl;
     cout << " W(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << length[e]<< endl;
     
     if (dual){
     abndncError += pow(double(abundance_map[g.target(e)] - abundance_map[g.source(e)]),2) / 10000;
     methylCallError += length[e];
     match++;
     }
     else{
     abndncError += pow(double(abundance_map[g.target(e)]),2) /10000 + pow(double(abundance_map[g.source(e)]),2) /10000;
     methylCallError +=1;
     }
     }
     int TP = match;
     int FN = trueMethylData.size() - match;
     int FP = estimatedMethylData.size() - match;
     
     methylCallError = methylCallError/(match + FN +FP);
     
     //evalFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/eval.txt", std::ios_base::app);
     
     //evalFile << "abndncError " << "\t"  << "methylCall Error " << "\t" << "TP "<< "\t" << "FN " << "\t" << "FP " << endl;
     evalFile << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
     
     //std:cout << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
     */
    
    
    
    
    
}


