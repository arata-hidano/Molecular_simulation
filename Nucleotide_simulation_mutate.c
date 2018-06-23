/* 
Nucleotide random forward-mutation model.
Overall Steps:
	P1: Reading sequence data and write output XML files
	P2: Add mutations
      (1) Start with a single/multiple isolates.
      (2) Mutation occurs according to Gillespie algorithm.
      (3) Randomly choose one isolate.
      (4) Randomly choose one nucleotide base.
      (5) Randomly choose one string to replace the base (note the exisiting string should be excluded in this random selection)
      (6) Replace and keep on.
      (7) When infection happens, give exactly the same set of sequence to infected one.
    P3: Add transmission: copying sequence of parent to that of child.

*/

/* Part 1: Simulating a nucleotide sequence*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <malloc.h>

/*Declare functions*/
void read_sequence_data() ;
int createdata();
void JukesCantor();

/*File names*/
char SequenceData[] = "FJ785304_1_VP1.txt"; //Seq file for FMD
char OutputXML[] = "test.xml"; //Output XML file name

/*Declare variables*/
#define seq_length_FMD 633 //length of VP1, FMD

char FMD_seq_master[634]; //Master (initial) sequence for FMD VP1
struct isolate {
	char sequence[634];
	int id;
	int parent_id ; //id of parent
	int d_time ; //day this isolate diverged from the parent
	int num_mut ; // number of mutations
};

/*=====================================================================================================================================================================*/
/* Start of main function*/

main(){


read_sequence_data(SequenceData,seq_length_FMD,FMD_seq_master) ;
printf ("<sequence taxon> %s </sequence>\n",FMD_seq_master);
system("pause");

struct isolate *master_isolate ;
master_isolate = (struct isolate*)malloc(sizeof( struct isolate )); 
master_isolate->d_time = 0;
master_isolate->id = 0;
master_isolate->num_mut = 0;
master_isolate->parent_id = 0;
strcpy(master_isolate->sequence, FMD_seq_master);
printf ("copied seq is %s", master_isolate->sequence);

system("pause") ;

/* Project 2: Mutations*/
JukesCantor(master_isolate,seq_length_FMD) ;
system("pause");

/* Project 3: Transmission*/

struct isolate *child_isolate ;
child_isolate = (struct isolate*)malloc(sizeof( struct isolate ));
copy_isolate(master_isolate,child_isolate,5,0) ;
printf("Copied isolate is %s", child_isolate->sequence);
printf("Copied isolate id is %d", child_isolate->id);
system("pause");
/// now copying isolate is DONE
createdata(OutputXML,FMD_seq_master,seq_length_FMD);

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
    seq_master[i] = '\0' ; // Note that in C string variable is created by " " and just character is by ''
    //system("pause");
              } 
}

/* Writing XML file*/


int createdata(char OutputXML[],char sequence[], int seq_length)
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
fclose(fb) ;
return 0 ;
}

/* Jukes-Cantor Mutations. Choose random position, find a replacing base, and replace*/

void JukesCantor(struct isolate *mutating_isolate, int seq_length)
{
	srand((unsigned)time(NULL)); //set seed
	int r_position = rand()%seq_length ; //random position at nucleotide seq (random number generator needs improvement)
	char next_base; // make sure it is just char, not char[], which becomes character string type
	char current_base ; //get the current base
    current_base = mutating_isolate->sequence[r_position] ;
//	strcpy(current_base, sequence[r_position]); //does not work because strcpy copies a string of characters
	printf("Current base is %c", current_base);
	
//	if(strncmp(current_base,'A',1) )
//	printf("Yes");
	
	if(current_base=='A') // make sure the character is bracketted with '', not " "
	{
		printf("Matched!");
		next_base= (char)"TGC"[rand()%3];
	
	}
	else if(current_base=='T')
	{
		printf("Matched!");
		next_base= (char)"AGC"[rand()%3];
		//strcpy(next_base, "AGC"[rand()%3]);
	
	}
	else if(current_base=='G')
	{
		printf("Matched!");
		next_base= (char)"ATC"[rand()%3];
		//strcpy(next_base, "ATC"[rand()%3]);
	
	}
	else if(current_base=='C')
	{
		printf("Matched!");
		next_base= (char)"ATG"[rand()%3];
		//strcpy(next_base, "ATG"[rand()%3]);
	
	}
	else
	{
		printf("Error! This is not eligible base!") ;
	}
	printf("Next base is %c", next_base);
	mutating_isolate->sequence[r_position] = next_base;
	printf("Now base is %c", mutating_isolate->sequence[r_position]);
}


/* Transmission of isolate*/
void copy_isolate(struct isolate *parent, struct isolate *child, int today, int current_id)
{
child->d_time = today;
child->id = current_id+1;
child->num_mut = 0;
child->parent_id = parent->id;
strcpy(child->sequence, parent->sequence) ; 
// we can add mutation that accompanies at mutation if needed	
	
}
