#include <iostream>
#include <fstream>
#include <string>
using namespace std;

string* joinfile(string sarray[]) // put each contig of the genome into an index of sarray
{
  string line;
  string linelist[905752]; // array is the length of lines in ham_G15
  ifstream myfile ("ham_pb4_shred_2k.utg.fasta");
  int counter = 0;
  if (myfile.is_open())
    {
      while ( getline(myfile,line))
	{
	      linelist[counter] = line;
	      counter++;
	}
      myfile.close();
    }
  int linecount = 0;
  int count = 0;
  while (linecount < 905752)
  {
      string sequenstr = ""; 
      if (linelist[linecount][0] == '>')
	{
	  sequenstr = linelist[linecount];
	  linecount = linecount + 1;
	  while ((linecount < 905752) && (linelist[linecount][0] !='>'))
	    {
	      sequenstr += linelist[linecount];
	      linecount = linecount + 1;
	    }
	  sarray[count] = sequenstr;
	  count++;
	}
      else
	{
	  linecount = linecount + 1;
	}
  }
  return sarray;
}

string* Ncontigs(string sarray[], string Ncontiglist[])
{
  ofstream myfile ("ham_pb4_shred_contig_N.txt");
  for(int i = 0; i < 6886; i++) // go through each contig
    {
      int begin = 0;
      string letter = "";
      int x = 0;
      for(int l = 0; l < sarray[i].size(); l++)
	{
	  if ((sarray[i][l] == '.') && (sarray[i][l+1] == '0') && (sarray[i][l+2] == '0'))
	    {
	      begin = l+3; // set begin to the start of DNA
	    }
	}
      for(int j = begin; j < sarray[i].size(); j++)
	{
	  if(sarray[i][j] == 'n' or sarray[i][j] == 'N') // if there is a "N" instead of a nucletoide in the DNA file
	    {
	      int counter = begin;
	      while (counter < sarray[i].size() - begin) // look at all the DNA, not including name of contig
		{
		  if (counter == begin) // first write name of contig, then first row of DNA
		    {
		      if (myfile.is_open())
			{
			  myfile <<"\n"+ sarray[i].substr(0,begin) + "\n"+ sarray[i].substr(counter,70);
			}
		    }
		  else
		    {
		      if (myfile.is_open())
			{
			  myfile <<"\n"+sarray[i].substr(counter,70);
			}
		    }
		  counter = counter + 70;
		}
	      int remainder = sarray[i].size() - counter;
	      if (remainder > 0)
		{
		if (myfile.is_open())
		  {
		    myfile << "\n"+sarray[i].substr(sarray[i].size()-begin-remainder, remainder);
		  }
		}
	      break;
	    }
	}
    }
  myfile.close();
  int sline = 0;
  string line;
  ifstream file ("ham_pb4_shred_contig_N.txt");
  if (file.is_open())
    {
      while (getline(file,line))
	{
	  Ncontiglist[sline] = line;
	  sline ++;
	}
    }

  file.close();
  return Ncontiglist;
}

string* NcontigList(string Nlist[], string Nconlist[])
{
  int num = 0;
    for(int x = 0; x < 78916; x++)
    {
      if (Nlist[x][0] == '>')
	{
	  Nconlist[num] = Nlist[x];
	  num++;
	}
    }

  return Nconlist;
}

bool compList(string ind, string sList[], int size) // compares the names of contigs
{
  for(int m = 0; m < size; m++)
    {
      if (ind == sList[m])
	{
	  return true;
	}
    }
  return false;
}

void NLessContigs(string Ncontigtitle[]) // get all contigs that do not have Ns
{
  string fullList[905752];
  int count = 0;
  int countNless = 0;
  string line;
  ifstream myfile("ham_pb4_shred_2k.utg.fasta");
  if (myfile.is_open())
    {
      while(getline(myfile,line))
	{
	  fullList[count] = line;
	  count++;
	}
      myfile.close();
    }
  int hamLen = 0;
  ofstream file("c++_ham_pb4_No_Ns.fasta");
  string myline;
  while (hamLen < 905752)
    {
      if (fullList[hamLen][0] == '>')
	{
	  if (compList(fullList[hamLen], Ncontigtitle, 708) == false) // if contig does not have any Ns write it to file
	    {
	      if (file.is_open())
		{
		  file << "\n" + fullList[hamLen];
		}
	      hamLen++;
	      countNless++;
	      while((hamLen < 905752) && (fullList[hamLen][0] != '>'))
		{
		  if (file.is_open())
		    {
			  file << "\n" + fullList[hamLen];
		    }
		  countNless++;
		  hamLen++;
		}
	    }
	  else
	    {
	      hamLen++;
	    }
	}
      else
	{
	  hamLen++;
	}
    }
}
int main()
{
  string sarray[6886];
  string Ncontiglist[156047];
  string Nconlist[708];
  string *strarray = joinfile(sarray); // put all contigs into a string array
  string *Nlist = Ncontigs(strarray,Ncontiglist); // get all the contigs that have "n" in them
  string *Ncontigtitle = NcontigList(Nlist, Nconlist); // get the names of the contigs with Ns in them
  NLessContigs(Ncontigtitle); // put contigs with no Ns in file
}
