/* 
Nucleotide random forward-mutation model.
Overall Steps:
      (1) Start with a single/multiple isolates.
      (2) Mutation occurs according to Gillespie algorithm.
      (3) Randomly choose one isolate.
      (4) Randomly choose one nucleotide base.
      (5) Randomly choose one string to replace the base (note the exisiting string should be excluded in this random selection)
      (6) Replace and keep on.
      (7) When infection happens, give exactly the same set of sequence to infected one.

*/

/* Part 1: Simulating a nucleotide sequence*/
#include <stdio.h>
#include <ctype.h>

/*Declare functions*/
void read_sequence_data() ;
void createdata();

/*File names*/
char SequenceData[] = "FJ785304_1_VP1.txt"; //Seq file for FMD
char OutputXML[] = "test.xml"; //Output XML file name

/*Declare variables*/
int seq_length_FMD = 633; //length of VP1, FMD
char FMD_seq_master[634]; //Master (initial) sequence for FMD VP1



/*=====================================================================================================================================================================*/
/* Start of main function*/

main(){


read_sequence_data(SequenceData,seq_length_FMD,FMD_seq_master) ;
printf ("<sequence taxon> %s </sequence>\n",FMD_seq_master);
system("pause");



createdata(OutputXML,FMD_seq_master,seq_length_FMD);
return 0;
}

/*END of main function*/

/*====================================================================================================================================================================*/
/* Functions*/

/* Reading sequence data*/

void read_sequence_data(char SequenceData[], int seq_length, char seq_master[])
{
     FILE *Seq = fopen(SequenceData, "r");
     int c;
     int i = 0;
     if (Seq) {
      while ((c = fgetc(Seq)) != EOF)
                {
    
                printf("c is %c\n", c);
                if (isspace(c)) // this was necessary because for some reasons the sequence text file contains empty character
                   {
                   printf("this is empty");
                   }
                else{
                     seq_master[i]=c;
                     printf("seq is %c\n", seq_master[i]);
                     //printf("i is %d\n", i);
                     i++;
                     }
                }
    fclose(Seq);
    //printf("i is %d\n", i);
    seq_master[i] = "\0" ;
    return(seq_master) ;
    //system("pause");
              } 
}

/* Writing XML file*/


void createdata(char OutputXML[],char sequence[], int seq_length)
{
     int i;
     FILE *fb=fopen(OutputXML,"w");
fprintf ( fb,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
fprintf ( fb,"<data dataType=>\n");
fprintf ( fb,"<sequence taxon>") ;
for (i=0;i<seq_length; i++){
	fprintf ( fb,"%c",sequence[i]);
	}
fprintf ( fb,"</sequence>\n") ;
fprintf ( fb,"</data>\n");
}

