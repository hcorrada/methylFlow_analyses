// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <algorithm>
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <lemon/list_graph.h>
//#include <lemon/lp_base.h>
#include <lemon/lp.h>
#include <lemon/maps.h>
#include <MethylRead.hpp>
//#include <lemon/bipartite_matching.h>
//#include <lemon/concepts/bpgraph.h>
//lemon::Lp lp;

using namespace methylFlow;
using namespace std;

int chr, var;


struct PosChr
{
    int pos;
    int chr;
    
    PosChr(int _p, int _c) : pos(_p), chr(_c) {}
    
    bool operator<(const PosChr& rhs) const
    {
        return ((chr < rhs.chr) || (chr == rhs.chr && pos < rhs.pos));
    }
    
    bool operator==(const PosChr& rhs) const
    {
        return (pos == rhs.pos && chr == rhs.chr);
    }
};

std::ifstream estimatedPatternFile, mFile, shortReadFile, samReadFile;
std::ofstream weightFile, matchFile, matchAppFile, mEstimatedPercentageFile, mReadPercentageFile, mSamPercentageFile, patternPruneSamFile, patternPruneReadFile ;
std::vector<MethylRead*> estimatedMethylData, readMethylData, samMethylData;
std::vector<float> estimatedAbundanceData, readAbundanceData, samAbundanceData;
std::vector<int> estimatedChrData, readChrData, samChrData, estimatedCidData;


std::set<PosChr> cpgPos;

std::map<PosChr, float> estimated_U_map, estimated_M_map, read_U_map, read_M_map, sam_M_map, sam_U_map;

std::map<int, int> minCidPos, maxCidPos;

void readEstimatedPattern(int start, int end){
    cerr << "reading patterns.tsv file" << endl;
    
	// "s" is for start position, "e" is for end position
    float abundance;
	string  methylString;
    int chr;
    string dummyLine;
    getline(estimatedPatternFile, dummyLine);
	int i = 0;
    int lastCid = -1;
	while(!estimatedPatternFile.eof()) {
        int  s = -1, e = -1, cid, pid;
        
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
            estimatedChrData.push_back(chr);
            estimatedCidData.push_back(cid);
            if (cid != lastCid) {
                lastCid = cid;
                minCidPos[cid] = s;
                maxCidPos[cid] = e;
                
            } else {
                minCidPos[cid] = min(minCidPos[cid],s);
                maxCidPos[cid] = max(maxCidPos[cid], e);
            }
            
		}
	}
    cout << "size of estimatedMethyl = " << i << endl;
    /*  float sum =0;
     int lastChr = 0;
     for(unsigned int j= 0; j<estimatedAbundanceData.size(); j++){
     if (estimatedChrData[j] == lastChr) {
     sum += estimatedAbundanceData.at(j);
     }
     else{
     for(unsigned int i= 0; i<j; i++){
     estimatedAbundanceData.at(i) = estimatedAbundanceData.at(i)/sum ;
     }
     lastChr = estimatedChrData[j];
     sum = estimatedAbundanceData.at(j);
     
     }
     
     }*/
    
    
    
}

void readShortRead(int start, int end){
    
    
    // need to define where we  get chr data from
    int  s, l;
	string readId, methylString, strand, etc;
    int chr = 1;
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
            readAbundanceData.push_back(1);
            readChrData.push_back(chr);
		}
	}
    
    /*float sum =0;
     int lastChr = 0;
     for(unsigned int j= 0; j<readAbundanceData.size(); j++){
     if (readChrData[j] == lastChr) {
     sum += 1;
     }
     else{
     for(unsigned int i= 0; i<j; i++){
     readAbundanceData.at(i) = readAbundanceData.at(i)/sum ;
     }
     lastChr = readChrData[j];
     sum = 1;
     
     }
     
     }*/
    
    
}

void readSamInput(int start, int end){
    const int READ_LIMIT = 100000000000;
    std::string input;
    
    cerr << "reading sam input" << endl;
    int count = 0;
    while (std::getline(samReadFile, input)){
        //std::cerr << "[methylFlow] 0 Discarding lines start with " << input << std::endl;
        if(input.size() && input[0] !='@') break;
        
    }
    while (count < READ_LIMIT) {
        
        if(count >0 )
            std::getline( samReadFile, input );
        count++;
        
        std::string readid, rStrand, methStr, substStr;
        std::string QNAME, RNAME, CIGAR, RNEXT, SEQ, QUAL, NM, XX, XM, XR, XG;
        int FLAG, POS, MAPQ, PNEXT, TLEN;
        int rPos, rLen;
        std::istringstream buffer(input);
        
        
        //parse SAM format
        buffer >> QNAME >> FLAG >> RNAME >> POS >> MAPQ >> CIGAR >> RNEXT >> PNEXT >> TLEN >> SEQ >> QUAL >> NM >> XX >> XM >> XR >> XG;
        if ( !buffer || !buffer.eof() ) {
            std::cerr << "[samEvaluation] Error parsing SAM input" << std::endl;
            break;
        }
        //parse chr name
        std::string chromosome;
        if (RNAME[0] == 'c' || RNAME[0] == 'C')
            chromosome =  RNAME.substr(3);
        else
            chromosome = RNAME.substr(0);
        chr = atoi(chromosome.c_str());
        
        MethylRead * m;
        
        rPos = POS;
        rLen = SEQ.length();
        //rLen = findLength(QNAME);
        //std::cout << "rLen " << rLen << std::endl;
        //  std::cout << "count " << count << std::endl;
        
        // construct object with read info
        m = new MethylRead(rPos, rLen);
        if (methStr != "*"){
            m->parseXMtag(XM, CIGAR);
            samMethylData.push_back(m);
            samAbundanceData.push_back(1);
            samChrData.push_back(chr);
        }
    }
    
    /*float sum =0;
     int lastChr = 0;
     for(unsigned int j= 0; j<samMethylData.size(); j++){
     if (samChrData[j] == lastChr) {
     sum += 1;
     }
     else{
     for(unsigned int i= 0; i<j; i++){
     samAbundanceData.at(i) = samAbundanceData.at(i)/sum ;
     }
     lastChr = samChrData[j];
     sum = 1;
     
     }
     
     }*/
    
    //sort(samMethylData.begin(), samMethylData.end(), CompareReadStarts());
    cerr << "reading sam input finished" << endl;
    
}


void computeMethylPercentage(){
    //cerr << "compute methyl percentage" << endl;
    for (unsigned int i=0; i < samMethylData.size(); i++) {
        for (unsigned int j=0; j < samMethylData.at(i)->cpgs.size(); j++) {
            int position = samMethylData.at(i)->cpgs.at(j).offset + samMethylData.at(i)->start();
            bool meth = samMethylData.at(i)->cpgs.at(j).methyl;
            int chr = samChrData.at(i);
            std::set<PosChr>::iterator it;
            
            if ((it = cpgPos.find(PosChr(position, chr))) == cpgPos.end()) {
                cpgPos.insert(PosChr(position, chr));
                sam_U_map[PosChr(position, chr)] = 0;
                sam_M_map[PosChr(position, chr)] = 0;
                estimated_U_map[PosChr(position, chr)] = 0;
                estimated_M_map[PosChr(position, chr)] = 0;
            }
            if (meth) {
                sam_M_map[PosChr(position, chr)] += 1;
                // cerr << samAbundanceData.at(i) << "\t" << "sam_M_map[position]= " << it->pos << "\t"<< sam_M_map[PosChr(position, chr)] << endl;
                
            }
            else{
                sam_U_map[PosChr(position, chr)] += 1;
                //cerr << samAbundanceData.at(i) << "\t" << "sam_U_map[position]= " <<  it->pos << "\t" << sam_U_map[PosChr(position, chr)] << endl;
                
            }
        }
        
    }
    
    // cout << "size of estimatedMethylData = " << estimatedMethylData.size() << endl;
    for (unsigned int i=0; i < estimatedMethylData.size(); i++) {
        //cerr << "estimated " << i << "\t" << estimatedMethylData.at(i)->cpgOffset.size() << endl;
        for (unsigned int j=0; j < estimatedMethylData.at(i)->cpgs.size(); j++) {
            int position = estimatedMethylData.at(i)->cpgs.at(j).offset + estimatedMethylData.at(i)->start();
            bool meth = estimatedMethylData.at(i)->cpgs.at(j).methyl;
            int chr = estimatedChrData.at(i);
            std::set<PosChr>::iterator it;
            
            if ((it = cpgPos.find(PosChr(position, chr))) == cpgPos.end()) {
                cpgPos.insert(PosChr(position, chr));
                sam_U_map[PosChr(position, chr)] = 0;
                sam_M_map[PosChr(position, chr)] = 0;
                estimated_U_map[PosChr(position, chr)] = 0;
                estimated_M_map[PosChr(position, chr)] = 0;
            }
            if (meth) {
                estimated_M_map[PosChr(position, chr)] += estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" << "estimated_M_map[position]= " <<  it->pos << "\t" << it->chr <<"\t"<< estimated_M_map[PosChr(position, chr)] << endl;
            }
            else{
                estimated_U_map[PosChr(position, chr)] += estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" <<  "estimated_U_map[position]= " << it->pos << "\t" << it->chr <<"\t" << estimated_U_map[PosChr(position, chr)] << endl;
                
            }
        }
        
    }
    
}

void writeMethylPercentageMatrix(){
    
    //sort(cpgPos.begin(), cpgPos.end());
    
    mEstimatedPercentageFile << "chr" << "\t" << "pos" << "\t" << "unmethyl" << "\t" << " methyl" << endl;
    std::set<PosChr>::iterator it;
    for (it=cpgPos.begin(); it!=cpgPos.end(); ++it){
        float sum = estimated_U_map[*it] + estimated_M_map[*it];
        if (sum > 0.0001)
            mEstimatedPercentageFile << it->chr << "\t" << it->pos << "\t"  << estimated_U_map[*it]/sum << "\t"  << estimated_M_map[*it]/sum << endl;
        else
            mEstimatedPercentageFile << it->chr << "\t" << it->pos << "\t"  << 0 << "\t"  << 0 << endl;
        
    }
    
    
    mSamPercentageFile << "chr" << "\t" << "pos" << "\t" << "unmethyl" << "\t" << " methyl" << endl;
    
    for (it=cpgPos.begin(); it!=cpgPos.end(); ++it){
        float sum = sam_U_map[*it] + sam_M_map[*it];
        if (sum > 0.0001)
            mSamPercentageFile << it->chr << "\t" << it->pos << "\t"  << sam_U_map[*it]/sum << "\t"  << sam_M_map[*it]/sum << endl;
        else
            mSamPercentageFile << it->chr << "\t" << it->pos << "\t"  << 0 << "\t"  << 0 << endl;
        
        
    }
    
}

float computeMethylPercentageError(set<PosChr> &pos, map<PosChr, float> &m1, map<PosChr, float> &u1, map<PosChr, float> &m2, map<PosChr, float> &u2, int minPos, int maxPos){
    float error = 0;
    std::set<PosChr>::iterator it;
    int n = 0;
    for (it=pos.begin(); it!=pos.end(); ++it){
        if (it->pos < minPos || it->pos > maxPos)
            continue;
        n++;
        float mm1 = m1[*it];
        float uu1 = u1[*it];
        float mm2 = m2[*it];
        float uu2 = u2[*it];
        
        if (mm1 ==0 && uu1==0) {
            if (mm2==0 && uu2 ==0) {
                error +=0;
            }
            else{
                error += pow(mm2/(mm2+uu2), 2);
            }
        }
        else{
            if (mm2 ==0 && uu2==0) {
                error += pow(mm1/(mm1+uu1),2);
            }
            else{
                error += pow(mm1/(mm1+uu1)- mm2/(mm2+uu2),2);
                
            }
        }
    }
    if (n == 0)
        return 0;
    else
        return sqrt(error/float(n));
}


void pruneEstimatedPattern(){
    bool check = false;
    vector<float> errorVec;
    patternPruneSamFile << "componentID"<< "\t"<<"PatternID" << "\t" << "methylPercentageError" << endl;
    
    //sort(estimatedAbundanceData.begin(), estimatedAbundanceData.end());
    int lastCid = -1;
    float error0;
    for (unsigned int i = estimatedAbundanceData.size() -1 ; i >= 0 ; i--) {
        int cid = estimatedCidData.at(i);
        if (cid != lastCid || i == 0) {
            if (lastCid != -1) {
                for (unsigned int i = 0; i < errorVec.size(); i++) {
                    patternPruneSamFile << lastCid << "\t" << i << "\t" << errorVec.at(i) << endl;
                }
                errorVec.clear();
            }
            error0 = computeMethylPercentageError(cpgPos, estimated_M_map, estimated_U_map, sam_M_map, sam_U_map, minCidPos[cid], maxCidPos[cid]);
            errorVec.push_back(error0);
            lastCid = cid;
            
        }
        
        for (unsigned int j=0; j < estimatedMethylData.at(i)->cpgs.size(); j++) {
            int position = estimatedMethylData.at(i)->cpgs.at(j).offset + estimatedMethylData.at(i)->start();
            bool meth = estimatedMethylData.at(i)->cpgs.at(j).methyl;
            int chr = estimatedChrData.at(i);
            if (meth) {
                estimated_M_map[PosChr(position, chr)] -= estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" << "estimated_M_map[position]= " <<  it->pos << "\t" << it->chr <<"\t"<< estimated_M_map[PosChr(position, chr)] << endl;
            }
            else{
                estimated_U_map[PosChr(position, chr)] -= estimatedAbundanceData.at(i);
                //cerr << estimatedAbundanceData.at(i) << "\t" <<  "estimated_U_map[position]= " << it->pos << "\t" << it->chr <<"\t" << estimated_U_map[PosChr(position, chr)] << endl;
                
            }
        }
        float error = computeMethylPercentageError(cpgPos, estimated_M_map, estimated_U_map, sam_M_map, sam_U_map, minCidPos[cid], maxCidPos[cid]);
        errorVec.push_back(error);
        if (error > 1.1 * error0 && check == false) {
            cout << "best Patterns are from 0 to "<<  i << endl;
            check = true;
        }
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
		start = atoi(argv[4]);
		end = atoi(argv[5]);
		var = atof(argv[6]);
        
	}
    
    std::string outdirname = argv[3];
    
    
    
    if (var == 0){
        cout << "your input is not a sam file" << endl;
        
        shortReadFile.open( argv[1]);
        buffer.str("");
        buffer << outdirname << "/methylPercentageRead.txt";
        mReadPercentageFile.open( buffer.str().c_str() );
        buffer.str("");
        buffer << outdirname << "/patternPruneRead.txt";
        patternPruneReadFile.open( buffer.str().c_str() );
        
    }
    else{
        cout << "your input is a sam file" << endl;
        
        samReadFile.open( argv[1] );
        buffer.str("");
        buffer << outdirname << "/methylPercentageSam.txt";
        mSamPercentageFile.open( buffer.str().c_str() );
        buffer.str("");
        buffer << outdirname << "/patternPruneSam.txt";
        patternPruneSamFile.open( buffer.str().c_str() );
    }
    
    
    estimatedPatternFile.open( argv[2] );
    
    buffer.str("");
    buffer << outdirname << "/methylPercentageEstimated.txt";
    mEstimatedPercentageFile.open( buffer.str().c_str() );
    
    
    
	//truePatternFile.open(argv[1]);
    //estimatedPatternFile.open(argv[2]);
    //evalFile.open(argv[3],std::ios_base::app);
    //MCFFile.open(argv[4],std::ios_base::app);
    
	
	//################     read the input .tsv data to the "line" number
	cout << "reading patterns " << start << endl;
    readEstimatedPattern(start , end);
    cout << "reading short read " <<  endl;
    
    if (var == 0) {
        readShortRead(start, end);
        computeMethylPercentage();
        cout << "write methyl Percentage  shortRead " << endl;
        writeMethylPercentageMatrix();
        mReadPercentageFile.close();
        patternPruneReadFile.close();
        
    }
    else{
        readSamInput(start, end);
        computeMethylPercentage();
        cout << "write methyl Percentage Sam " << endl;
        writeMethylPercentageMatrix();
        mSamPercentageFile.close();
        patternPruneSamFile.close();
    }
    
    
    
    
    
}


