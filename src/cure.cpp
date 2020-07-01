//Merge multiple molecule dataset by clean and fusion all datasets.
//Jintao Meng

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <set>
#include <map>
#include <queue> 
#include <deque>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

typedef pair<int, float> pif; 
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef vector<string> vs;
typedef vector<vs> vvs;
typedef vector<float> vf;
typedef map<string, vf> mvf;

typedef pair<float, float> pff;
typedef vector<pff> vp;
typedef map<string, vp> mvp;

#define debug(x)	cout << #x << " = " << x << endl
#define sz			size()
#define mod(A,B)	((((A) % (B)) + (B)) % (B))
#define b2e(A)		(A).begin(), (A).end()
#define e2b(A)		(A).rbegin(), (A).rend()
#define rev(A)		 reverse(b2e(A))
#define s(A)		sort(b2e(A))
#define ss(A)		stable_sort(b2e(A))
#define Fill(A,B)	fill(b2e(A), B)
#define minelt(A)	*min_element(b2e(A))
#define maxelt(A)	*max_element(b2e(A))
#define contains(A,B)   (find(b2e(A),B)!=(A).end()) 

#define un(A)		unique(b2e(A))
#define er(A)		(A).erase(un(A), (A).end())

#define REP(i,n)	for(int i=0;i<(n);++i)
#define EACH(i,x)	REP(i,(x).size())
#define FOR(i,a,b)	for(int i=(a);i<(b);++i)

#define Min(a,b)	((a)<(b)?(a):(b))
#define Max(a,b)	((a)>(b)?(a):(b))

vs tok(string a){vs r; char *b=new char [a.length()+1],*s; strcpy(b,a.c_str()); for(s=strtok(b," ");s;s=strtok(0," ")) r.push_back(s); return r;}
vs tokS(string a, string d){vs r; char *b=new char [a.length()+1],*s; strcpy(b,a.c_str()); for(s=strtok(b,d.c_str());s;s=strtok(0,d.c_str())) r.push_back(s); return r;}

void printMolData(string name, mvp data)
{
	cout<<"Dataset: "<<name<<endl;
	for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{
		cout<<it->first;
		REP(j,it->second.size())
			cout<<','<<it->second[j].first<<'-'<<it->second[j].second;
		cout<<endl;
	}   
}

void writeMolData(string filename, mvp data, int orgflag=0)
{
	ofstream out(filename.c_str());

	if (!out.is_open())
	{
		cout<<"Failed to open "<<filename<<endl;    
		return;
	}    
	
	if(orgflag==0)		out<<"smiles,logS,weight"<<endl;
	else			out<<"smiles,logS"<<endl;
	for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{
		REP(j,it->second.size())
			if(orgflag==0)	out<<it->first<<','<<it->second[j].first<<','<<it->second[j].second<<endl;
			else		out<<it->first<<','<<it->second[j].first<<endl;
	}   
	out.close();
}

pif calcMultiMolData(mvp data, float delta=0)
{
    	int cnt=0;
    	float percent;
	cout<<"------------"<<endl;
	for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{
		int times = it->second.size();
		int max=-1000000, min=10000000;
		for(vp::iterator vt=it->second.begin(); vt!=it->second.end(); vt++)
		{	
			if(vt->first>max)	max=vt->first;
			if(vt->first<min)	min=vt->first;
		}
	        float diff = max-min;  	
		if(diff<delta)	
		{	
			times = 1;
		}
		else 
		{		
			cout<<it->first;		
			REP(j,it->second.size())
				cout<<','<<it->second[j].first<<'-'<<it->second[j].second;
			cout<<endl;
		}
		if(times>1)	cnt++;
	}	
	percent = 100.0*cnt/data.size();
	return pair<int, float> (cnt, percent);
}

mvp cleanParseData(mvp data, int norm=0, float dataWeight=1.0)
{
        string split = string(".,U,Ge,Pr,La,Dy,Ti,Zr,Rh,Lu,Mo,Sm,Sb,Nd,Gd,Cd,Ce,In,Pt,Sb,As,Ir,Ba,Hg,Se,Sn,Ti,Fe,Si,Al,Bi,Pb,Pd,Ag,Au,Cu,Pt,Co,Ni,Ru,Mg,Zn,Mn,Cr,Ca,K,Li,SF5,SF6");
        vector<string> rsplit = tokS(split, string(","));

        mvp ret;
        for(mvp::iterator it=data.begin(); it!=data.end(); it++)
        {
                int i, flag = 0;
                for(i=0;i<rsplit.size();i++)
                        if(it->first.find(rsplit[i])!=string::npos)     {flag=1; break;}
                if(flag)
                {
                	cout<<"Containing: "<<rsplit[i]<<endl;
                        cout<<"Removing: "<<it->first<<endl;
                        continue;
                }
                ret[it->first] = it->second;
        }
        for(mvp::iterator it=ret.begin(); it!=ret.end(); it++)
                sort(it->second.begin(), it->second.end());

        for(mvp::iterator it=ret.begin(); it!=ret.end(); it++)
	{
		vp tmp(1, it->second[0]);	
		for(int i=1; i<it->second.size(); i++)
		{
			if(tmp[tmp.size()-1].first==it->second[i].first)
				tmp[tmp.size()-1].second += it->second[i].second;
			else	tmp.push_back(it->second[i]);				
		}
		it->second = tmp;
	}	
	
	if(norm)
        for(mvp::iterator it=ret.begin(); it!=ret.end(); it++)
        {
        	float wegt = 0;
		for(int i=0; i<it->second.size(); i++)
			wegt += it->second[i].second;
		for(int i=0; i<it->second.size(); i++)
			it->second[i].second /= wegt;	
		for(int i=0; i<it->second.size(); i++)
			it->second[i].second *= dataWeight;	
        }

	return ret;	
}

vp weightedVoteWithSlidingWindow(string smiles, vp value, float slide=0)
{
	if(value.size()==1)	return value;

        float max=-100000, min=100000;
        for(int head=0;head<value.size();head++)        //[end, head]
        {
                if(value[head].first>max)       max=value[head].first;
                if(value[head].first<min)       min=value[head].first;
        }
	int debug = 0;
	if(max-min>slide)	debug=1;	

	debug=0;
	if(debug)
	{
        cout<<smiles<<":";
	for(int head=0;head<value.size();head++)	//[end, head]
		cout<<value[head].first<<"#"<<value[head].second<<"  ";
	cout<<endl;
	}

	vp val = value;
	float dist = 0;
	while(1)
	{
		float localdist = 10000;
		int index;

		for(int i=0;i<val.size()-1;i++)
		{
			if(localdist > fabs(val[i].first-val[i+1].first))		
			{
				localdist = fabs(val[i].first-val[i+1].first);
				index = i;
			}
		}

		if(localdist > slide)	break;
			
		val[index].first = val[index].first*val[index].second + val[index+1].first*val[index+1].second;	
		val[index].second = val[index].second + val[index+1].second;
		val[index].first /= val[index].second;
		val.erase(val.begin() + index +1);
		
		if(0)
		{
		cout<<"after merge"<<endl;
		for(int i=0;i<val.size();i++)
			cout<<val[i].first<<"#"<<val[i].second<<" ";
		cout<<endl;	
		}
	}
	if(debug)
	{
//        cout<<smiles<<":";
	REP(j,val.size())
		cout<<val[j].first<<'#'<<val[j].second<<" ";
	cout<<endl;
	}
	return val;
}

mvp cureData(vector<mvp> datasets, mvp data, float slide=0.5, float threshold=1.0)
{
        mvp collection;
        REP(i,datasets.size())
        for(mvp::iterator it=datasets[i].begin(); it!=datasets[i].end(); it++)
        {
                string key = it->first;
                for(vp::iterator vt=it->second.begin(); vt!=it->second.end(); vt++)
                {
                        float logS = vt->first;
                        float wegt = vt->second;
                        if(collection.find(it->first)==collection.end())
                                collection[it->first] = vp(1,pff(logS,wegt));
                        else    collection[it->first].push_back(pff(logS, wegt));
                }
        }

        for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
                sort(it->second.begin(), it->second.end());

        for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
                it->second = weightedVoteWithSlidingWindow(it->first, it->second, slide);

        for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{	
//		it->second = vp(1, pff(0.0,0.0));
		it->second = collection[it->first];
	}

	//rewrite here. 
	float maxWeight=0;		
	for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{	
		float tot=0;
		REP(j,it->second.size())
			tot += it->second[j].second;
		if(tot>maxWeight)	maxWeight = tot;
	}

        cout<<"maxWeight"<<maxWeight<<endl;	
	if(maxWeight>threshold)	maxWeight = threshold;
	for(mvp::iterator it=data.begin(); it!=data.end(); it++)
	{
		float tot=0;
		REP(j,it->second.size())
			tot += it->second[j].second;
		
		REP(j,it->second.size())
		if(tot<maxWeight)	it->second[j].second /= maxWeight;
		else			it->second[j].second /= tot;
	}

        return data;	
}

mvp mergeMultiData(vector<mvp> datasets, float slide=0.5)
{
        mvp collection;	
    	REP(i,datasets.size())
	for(mvp::iterator it=datasets[i].begin(); it!=datasets[i].end(); it++)
	{
		string key = it->first;
		int    cnt = it->second.size();
                for(vp::iterator vt=it->second.begin(); vt!=it->second.end(); vt++)
		{
			float logS = vt->first;
                        float wegt = vt->second/cnt;
			if(collection.find(it->first)==collection.end())	
				collection[it->first] = vp(1,pff(logS,wegt)); 
			else	collection[it->first].push_back(pff(logS, wegt));							
		}		
	}

	for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
		sort(it->second.begin(), it->second.end());	
        
	for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
		it->second = weightedVoteWithSlidingWindow(it->first, it->second, slide);	

	float maxWeight=0;		
	for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
		if(it->second[0].second>maxWeight)	maxWeight = it->second[0].second;

	for(mvp::iterator it=collection.begin(); it!=collection.end(); it++)
		it->second[0].second /= maxWeight;

	return collection;	
}

pif calcMolDataAinB(mvp dataA, mvp dataB)
{
    	int cnt=0;
    	float percent;

	for(mvp::iterator it=dataA.begin(); it!=dataA.end(); it++)
		if(dataB.find(it->first)!=dataB.end())	cnt++;
	percent = 100.0*cnt/dataA.size();
	return pair<int, float> (cnt, percent);
}

pair<float, float> calcMolAB(mvp dataA, mvp dataB)
{
    	int   cnt_same=0, cnt_diff=0;
    	float pct_same=0, pct_diff=0;

	int   total = 0;	
	for(mvp::iterator ita=dataA.begin(); ita!=dataA.end(); ita++)
	{	
		total += ita->second.size();
		mvp::iterator itb = dataB.find(ita->first);
		if(itb==dataB.end())	continue;	
		for(int i=0; i<ita->second.size(); i++)
		{
			int dis_same = 0;
			int dis_diff = 0;
			for(int j=0;j<itb->second.size();j++)
			{
				float dis = fabs(ita->second[i].first - itb->second[j].first);
				if(dis<=0.01)	dis_same = 1;
				else		dis_diff = 1;		
			}
			cnt_same += dis_same;
			cnt_diff += dis_diff;
		}
	}
	pct_same = 100.0 * cnt_same/total;
	pct_diff = 100.0 * cnt_diff/total;
	return pair<float, float> (pct_same, pct_diff);
}


vector<mvp> ReadDatasets(vs files, vs names)
{
    vector<mvp> datasets;

    REP(i,files.size())
    {
        string line;
        ifstream myfile(files[i].c_str());

        if (!myfile.is_open())
        {
                cout<<files[i]<<" file is not exist"<<endl;
                continue;
        }

        mvp data;
        int cnt = 0;
        while ( getline (myfile,line) )
        {
                float weight=0;
                if(!(cnt++))    continue;
                int sepcnt = count(line.begin(), line.end(), ',');
                vs table = tokS(line, string(","));
                if(sepcnt==0)           continue;
                else
                {
                        string label = table.size()>2?table[2]:string("");
                        string temp  = table.size()>3?table[3]:string("");

                        if(names[i].find("aqua")!=names[i].npos)        weight = 0.9;
                        else if(names[i].find("esol")!=names[i].npos)   weight = 0.9;
                        else if(names[i].find("phys")!=names[i].npos)   weight = 1.0;
                        else if(names[i].find("moe")!=names[i].npos)    weight = 0.9;
                        else if(names[i].find("ochem")!=names[i].npos)
                        {
                                weight = 0.85;
                                if(label.find("Training")!=label.npos)  weight = weight * 1.0;
                                else if(label.find("Test")!=label.npos) weight = weight * 0.5;
                                else                                    weight = weight * 0.2;
                        }
                        else if(names[i].find("aqsol")!=names[i].npos)
                        {
                                weight = 0.6;
                                if(label.find("G1")!=label.npos)        weight *= 1.0;
                                else if(label.find("G2")!=label.npos)   weight *= 0.4;
                                else if(label.find("G3")!=label.npos)   weight *= 0.6;
                                else if(label.find("G4")!=label.npos)   weight *= 0.2;
                                else if(label.find("G5")!=label.npos)   weight *= 0.2;
                                else                                    weight  = 0;
                        }
                        else if(names[i].find("chembl")!=names[i].npos)
                        {
                                weight = 0.70;
                                if(label.find("pH 7")!=label.npos)      weight *= 1.0;
                                else if(label.find("pH 6")!=label.npos) weight *= 1.0;
                                else                                    weight *= 0.5;
                        }
                        else if(names[i].find("kinect")!=names[i].npos)
			{
				weight = 1.0;
				if(label==string("") || temp==string(""))	weight *= 0.5;
				else
				{
                                	float pHvalue = std::stof(label);
					float tmpvalue = std::stof(temp);
			 		if(pHvalue>6-0.001 && pHvalue<8+0.001)	weight *=1.0;
					else					weight *=0.5;
				
					if(tmpvalue>20-0.001 && tmpvalue<30+0.001)	weight *=1.0;
					else						weight *=0.5;			
				}
			}
                        else            weight = 0;
                }
                if(weight==0)   continue;
                if(data.find(table[0])==data.end())     data[table[0]] = vp(1, pff(atof(table[1].c_str()), weight));
                else                                    data[table[0]].push_back(pff(atof(table[1].c_str()), weight));
        }

        datasets.push_back(data);
        myfile.close();
    }

    return datasets;
}

int main(int argc,char *argv[])
{
    vector<mvp> datasets;
    vs files, names;
    FOR(i,1,argc)		files.push_back(string(argv[i]));
    FOR(i,0,files.size())	names.push_back(tokS(files[i], string("/")).back());

    datasets = ReadDatasets(files, names);

//    REP(i,files.size())
//	writeMolData(names[i]+string("org"), datasets[i], 1);					

    REP(i,files.size())
	writeMolData(names[i], cleanParseData(datasets[i], 1, 1.0));					

    REP(i,files.size())
    {
	float wedg=1;		
	if(names[i].find("aqua")!=names[i].npos)	wedg=0.9;
	else if(names[i].find("esol")!=names[i].npos)	wedg=0.9;
	else if(names[i].find("phys")!=names[i].npos)	wedg=1.0;
	else if(names[i].find("moe")!=names[i].npos)	wedg=0.9;
	else if(names[i].find("ochem")!=names[i].npos)	wedg=0.85;
	else if(names[i].find("aqsol")!=names[i].npos)	wedg=0.5;
	else if(names[i].find("chembl")!=names[i].npos)	wedg=0.7;	
	else if(names[i].find("kinect")!=names[i].npos)	wedg=1.0;	
        else 	wedg = 0;

	//todo: Normalization 
	datasets[i] = cleanParseData(datasets[i], 1, wedg);
    }

    //calculate percent of multiple molecules in each dataset
    REP(i,files.size())
    {
	cout<<files[i]<<endl;
    }

    //calculate percent of molecules from A contains in B
    cout<<"\t\t";
    REP(i,files.size())	cout<<names[i]<<"\t";
    cout<<endl;
    REP(i,files.size())
    { 	
	cout<<names[i]<<"\t";
	REP(j,files.size())
	{
//		pif ret = calcMolDataAinB(cleanParseData(datasets[i]), cleanParseData(datasets[j]));
		pair<float, float> ret = calcMolAB(cleanParseData(datasets[i]), cleanParseData(datasets[j]));
		cout<<ret.first<<'-'<<ret.second<<"\t";
	}
	cout<<endl;
    }

//    mvp all = mergeMultiData(datasets, 0.5);
    mvp all = cureData(datasets, datasets[datasets.size()-1], 0.5, 2);
    writeMolData(string("collection.csv"), all);					

    return 0;      
}

