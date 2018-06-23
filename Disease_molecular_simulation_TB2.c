/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
What's new?
7th Nov 2017: Adding slaughter 
9th Nov 2017: Adding TB disease = 0 codes
13th Nov 2017: Parameter set 0 TB: movement ban 30days, max_S3toS4_TB = 365, no recovery. 
14th Nov 2017: Parameter set 1 TB: movement ban 90days, max_S3toS4_TB = 365*2
14th Nov 2017: Parameter set 0 FMD: movement ban 30 days 
18th Nov 2017: Number of isolate limits between 50 and 150
19th Nov 2017: Only one sample per farm allowed. Change FarmDataStatus[i][5] setting.
30th Mar 2018: Add subsampling in each sampling scenario after the end of each iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/




/* Disease-molecular simulation models

(1) Demographic of populations
1. Each herd has 3 production types: calf, heifer, adult. This is set as an array.
2.  in each production type, struct of animals are connected.
3. These structs move between production types if necessary.
4. Struct also moves between farms.

=====================================================
== Data =============================================
There are two pointers of farm management unit.
1. Pointer to struct of animals.
2. Pointer to data storing infection status and other info of farms.


Farm data (dynamic memory)
Animal (dynamic memory)
Movement (dynamic memory)

(2) Disease transmission model
1. Disease spread between farms via movement of infected animals or within a farm via contact with infected animals
2. Need parameters for transmission coefficient for FMD and TB

(3) Sequence mutation model
1. Mutation occurs at a constant rate, which will be included in Gillespie algorithm

*/
/*======================================================================================
===========Definition, variable declaration============================================*/

/* C LIBRARIES TO INCLUDE */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <ctype.h>

/*DEFINE VARIABLES*/
#define num_production_type 3
#define num_total_farms 45965  // this is the total number of herds that appeared in movement\longvity table
#define num_eligible 12213 //toy2 3672// 159744 if use toy1 // this is the number of eligible animals for seeding provided those farms had out-movements of adults and these animals are adult
#define  n_rid 16
int n_rid2 = n_rid ;
int num_total_farms2 = num_total_farms ;
int const num_animals = 3624420;  // as of 2000 July 1st
int const num_total_animals = 16534951;  
#define sim_years 10 //@@@ 
int num_max_isolate = 100000 ;  // should change this
#define seq_length_FMD 633 
#define seq_length_TB 3500  //@@@length of TB gene to be checked
int seq_length ;
int *died_out;
int *sample_limit ;
/*VARIABLES TO CHANGE PARAMETERS*/

/*===========================================================================================================*/
/*Change parameter depepnding on the disease*/
int disease = 1; //0 if TB, 1 if FMD
struct isolate {
		char sequence[seq_length_FMD+1]; //change this accordingly
		int id; // id of the isolate
		int parent_id ; //id of parent
		double d_time ; //day this isolate diverged from the parent
		double isolated_time; //day when it's isolated
		int num_mut ; // number of mutations
		long long akey ;
		int isolated; //whether it was isolated 1 is yes
		long long pro_id ; //which production unit it came from
		
					};

/*MUTATION MODEL OR NOT?*/
int mutation_model = 1; //1 if mutation model
int detection_scenario = 0 ; //0 is unbiased 1 is biased sampling (area with identified case higher detection)
int num_isolate_upper = 150;
int num_isolate_lower = 50;
int tot_iterations = 30;
/*============================================================================================================*/

long long num_farm_production = num_total_farms*num_production_type ; // This is the ID for each farm production unit
int num_farm_var = 6; // 11 variables for each farm; [0] farm_id, [1] x, [2] y, [3] testarea_num, [4] disease_status, [5] testing schedule, [6] tb_detected, [7] sum of occult(reactive) for TB, infectious for FMD, [8] sum of infectious for TB, immune for FMD
  int num_farm_var2 = 6;                      // [9] farm id that infected this farm, [10] date that this source farm infected this destination farm
// [0]disease status [1] detected status [2] sum R3 [3]sum R4 [4]test_year
char FarmDataFile[] = "all_herds_details_export2.csv";
int column_prostatus = 9 ; // number of variable in the table of infection status
long long current_pro_id ;
// at this monemnt; [0] N, [1] sus, [2] exposed, [3] occult, [4] detectable,[5] infection pressure, [6] cumu_pressure, [7] mutation_pressure, [8] move_ban

/* Variables related to individual animal*/

int num_animals_var = 9; // [0] akey, [1] farm_id, [2] age_day, [3] age_type, [4] sex, [5] breed, [6] pregnant_status, [7] disease_status, [8] eligible for seed
// just need to initialise [7] and [8]
char AnimalDataFile[] = "tanimals_present_1july2000_toy3.csv"; //toy2 has initial infection for calf, toy for adults.

/* Variables related to events(movement, birth,test)*/
char MoveDataFile[] = "tmovements_since2000july_export3.csv";
//char MoveDataFile[] = "tmovements_since2000july_export2.csv"; ///need to make it shor
//long num_moves = 9681150;
//long num_moves = 8842666; // this is from 1JUly2000
//long num_moves = 8840836;  //this is from 2JULY2000
long num_moves = 8860567 ; //this is from 2JULY2000 AND WITH MISSING VALUES
int num_moves_vars = 6; // serial_akey, date, src_serial_hid, des_serial_hid, age_type, src_testarea


char BirthDataFile[] = "birth_table_since2000july.csv";
int num_births = 9667100 ;
int num_births_vars = 5;// akey, bdate, src_farm, sex, breed 

/* //@@@ IGNORE TEST
char TestDataFile[] = "/C_run/tb_test_schedule_export.csv";*/
//int num_tests = 742 ;
//int num_tests_vars = 2;
int random_year ;

//slaughter
char SlaughterDataFile[] = "serial_akey_culled.csv";
long num_slaughter = 6165850 ;
int num_slaughter_vars = 2 ;

//community membership
char MembershipDataFile[] = "community_membership_2000_2010.csv";
char hid_rid[] = "hid_rid.csv";
char num_active_herds_rid[] = "num_active_herds_rid.csv" ;
int total_active_herds = 18706 ;
/* Variables related to simulations*/

int sim_days= 365*sim_years; // X years
//341 farms detected TB in this 3 years
long long i, j;
int iteration = 0;
int num_OutPut = 2; //infected and detected
char OutPutFile[] = "OutPut2.csv";
/* Variables related to disease*/
int day_S3toS4, day_S2toS3;
int max_S3toS4_TB = 365*2; // max day to detect - 1 Conlan(2014) modified
int max_S2toS3_TB = 0.4*365; // max day to occult - 1  Brooks-Pollock's parameter 11.1
int min_S3toS4_TB = 60 ;

int max_S3toS4_FMD = 300 ;
int min_S3toS4_FMD = 60 ;
int max_S2toS3_FMD = 100;

int max_S2toS3, max_S3toS4, day_to_S3, day_to_S4,min_S3toS4;
int movement_ban_days = 180 ;
double Se_S3 = 0.5;
double Se_S4 = 0.85;
double Se_slaughter = 0.8 ;// @@@@ have to revise this value
double beta_a_FMD = 5.0;//within-herd transmission parameter
double beta_a_TB = 0.05 ;
double beta_a ;
double beta_b = 0;//wildlife transmission parameter
double mu_FMD = 0.012*633/365.0; //previous value 0.012 // value from Di Nardo 2014 disease==1 use substitution rate used in Hall, which was further converted to the substitution rate for the whole VP1 seq and per day
//7.4*10E-6 //2.7*10-3 was the original
double mu_TB = 0.001*seq_length_TB/365.0 ;
double mu ;
double detection_pressure_origin ;
double detection_pressure_origin1 = 0.3/30.0;
double detection_pressure_origin0 = 0 ;
double detection_pressure_increased = 0.6/30.0;
double detection_pressure_decreased = 0.1/30.0; 

int rid_status[n_rid] ;

int num_active_herds_rid_array[n_rid] = {0} ;
 //must be random draw from some distribution 0.003 - 0.036 Brooks-Pollock
int *current_id_isolate; // a pointer to a memory that keeps the track of isolate id


/*File names*/
char SequenceData[] = "FJ785304_1_VP1.txt"; //Seq file for FMD
char OutputXML[] = "test.xml"; //Output XML file name
const char* OutIsolateDataFile_TB = "TB_Sequence_data_iteration";
const char* OutNexusFile_TB = "TB_Nexus";
const char* OutIsolateDataFile_FMD = "FMD_Sequence_data_iteration";
const char* OutNexusFile_FMD = "FMD_Nexus";

//char* OutIsolateDataFile_origin ;//"Sequence_data_iteration";
const char* extension = ".csv" ; 
//char* OutNexusFile_origin ;// "Nexus";
const char* extension_nexus = ".nex" ;
const char* extension_nexus_E = "E.nex" ;
const char* extension_nexus_R = "R.nex" ;
const char* extension_nexus_H = "H.nex" ;
const char* extension_nexus_H_I = "H_I.nex" ;

int num_isolate_vars = 8;
int num_isolate_seq_vars = 2;
// [0]: akey, [1]: pro_id, [2]: isolated_time, [3]: diverged_time, [4]: id, [5]: parent_id, [6]: num_mut, [7]: isolated (1==Yes).
//int *nth_isolate;
int *nth_isolate_detected ;
int *nth_slaughter;
int *nth_between ;
int *max_isolate_reached ;

/*Counter for Markov event*/
int *c_transmission ;
int *c_mutation ;
int *c_detection ;

/* STRUCTURE DECLARATIONS */  
struct animal_node {
      long long akey ; // animal id
      int age_day;   /* age in months or day*/
      int type;         // production type, 0 = calf, 1 = heifer, 2 = adult
      int sex;          /* sex, female = 1 male = 0*/
      int breed;
      int disease_status;    /* disease status. Sus = 0, Exposed = 1, Occult =2, Detectable =3. Occult and Detectable are infectious*/
      int pregnant_status; /* Pregnant status positive = 1*/
      int num_births; // record how many births it gave to
      long long current_pro_id ;
      struct animal_node *previous_node ;/*pointer to the previous node*/
      struct animal_node *next_node; /* pointer to next animal*/
      struct isolate *ptr_isolate ; //pointer to isolate information
   }; 
   
struct event_node {
      int event_type;   
	  /*event type: 0 calf movement
	                1 heifer movement
					2 adult movement
					3 new birth
					4 cull/death
					5 test 
					6 exposed to occult
					7 occult to detectable
					8 Unit change
					9 Mutation*/
      long long akey ; //animal id
      long long src_pro_id;  // production unit id where the event occurs, if it's movement then it's source farm
      long long des_pro_id; 
	  int src_testarea;         
      struct event_node *next_node; /* pointer to next event*/
   }; 
   
   
struct inf_move {
	long long des_farm_id ;
	int date ;
	struct inf_move *next_node ;
};

	

/* FUNCTION DEFINITIONS */
  
void read_farm_data(); //Function for reading the data
void read_movement_data() ;
void read_animal_data();
void read_test_data();
void read_membership_data() ;
void add_event_node();
void add_animal_node();
void visualize_list();
void visualize_animals() ;
void move_animal_unit();
void read_birth_data();
void read_slaughter_data() ;
void read_hid_rid();
void read_num_active_herds_rid() ;

void test_farms();
double update_markov_date();
void count_farms_infected_detected();
int write_OutPut();
void move_production_type();
 
void read_sequence_data() ;
int createdata();
void JukesCantor();
int export_seq_nexus() ;
int export_IsolateData() ;
/* Variables related to farm and FarmProductionStatus*/
//int num_farms = 16950; // number of farms present at 2000 july 




/*-------START OF THE MAIN PROGRAMME------------------------------*/
/* MAIN PROGRAM */
int main(void){
char* OutIsolateDataFile_origin ;
char* OutNexusFile_origin ;
  if(disease==0) //TB
      	{
      		max_S2toS3 = max_S2toS3_TB;
      		max_S3toS4 = max_S3toS4_TB ;
      		min_S3toS4 = min_S3toS4_TB;
      		beta_a = beta_a_TB ;
      		seq_length = seq_length_TB + 1;
      		mu = mu_TB ;
      		OutIsolateDataFile_origin = malloc(strlen(OutIsolateDataFile_TB)+1) ;
      		OutNexusFile_origin = malloc(strlen(OutNexusFile_TB)+1) ;
      		strcpy(OutIsolateDataFile_origin, OutIsolateDataFile_TB);
			strcpy(OutNexusFile_origin, OutNexusFile_TB);
      		
      	}
		else if(disease==1)
		{
			max_S2toS3 = max_S2toS3_FMD ;
			max_S3toS4 = max_S3toS4_FMD ;
			min_S3toS4 = min_S3toS4_FMD;
			beta_a = beta_a_FMD ;
			seq_length = seq_length_FMD + 1 ;
			mu = mu_FMD ;
			OutIsolateDataFile_origin = malloc(strlen(OutIsolateDataFile_FMD)+1) ;
      		OutNexusFile_origin = malloc(strlen(OutNexusFile_FMD)+1) ;
			strcpy(OutIsolateDataFile_origin, OutIsolateDataFile_FMD);
			strcpy(OutNexusFile_origin, OutNexusFile_FMD);
		
		  } 

char Seq_master[seq_length]; 


       
    
/*=========================DATA PREPARATION=====================================================*/      
       
    /* 2. Read Data*/
	/* Following creates data or structure as follows.
	(1) FarmData (stores farm variable e.g. xy, bTB status, testing area, farmer type)
	(2) Farm_pointer (stores pointers to the first animal on a given farm with a given management type)
	(3) Read animal data
	(4) MoveData (Movement data)
	(5) TestData(testing and sampling data)
	(6) BirthData (new born animals)
	(7) sequence data
	*/  
/*==========2.1  Read in Farm Data=============================================================== */
    	// This is the basic farm data such as number of cattle in each management status
   	 	/* First create array FarmData, which specifies the row memory size*/
   	 	/*2.8 Read community_membership_data*/			
   	  	 double **FarmData =malloc( sizeof(double *) * num_total_farms2);
   	  	 int **MembershipData = malloc(sizeof(int*)*num_total_farms2) ;
   	  	 for(i = 0; i < num_total_farms; i++)
      	  {
     	       /*Then, for each row, specifies the required memory size for column*/
     	     FarmData[i] = malloc( sizeof(double) * num_farm_var); 
			 MembershipData[i] = malloc(sizeof(int)*(sim_years)) ;
     	   }
      	
     	read_farm_data(FarmDataFile, FarmData, num_total_farms2);
     	read_membership_data(MembershipDataFile,MembershipData,num_total_farms2) ;
		 //printf("farm data read");     
/*=====2.3 Read in animal Data===========================================================*/
       		long long **AnimalData = (long long**)malloc( sizeof(long long *) * num_animals);
   	  		 for(i = 0; i < num_animals; i++)
      	 	 {
     	       /*Then, for each row, specifies the required memory size for column*/
     	     AnimalData[i] = (long long*)malloc( sizeof(long long) * num_animals_var);      	   		
			}
					
      	read_animal_data(AnimalDataFile,AnimalData,num_animals) ;
      //	printf("The first animals is %lld",AnimalData[0][0]) ;
      	//printf("animal data read");
      /* Check if it reads AnimalData properly
	  	printf("Animal id is %lld",(long long) AnimalData[0][0]); // intersting it's now working, why do we need declaration?
        long long temp = AnimalData[0][0];
        int temp2 = AnimalData[0][1];
        printf("Animal id and farm is %lld, %d", temp, temp2);
      	system("pause") ;*/
//printf("b") ;
/*2.4 Read movement data*/
   	  	 long long **MoveData = (long long**)malloc( sizeof(long long *) * num_moves);
   	    	for(i = 0; i < num_moves; i++)
   	    	 {
   	    	 MoveData[i] = (long long*)malloc( sizeof(long long) * num_moves_vars);  
   		     	}
   	    	read_movement_data(MoveDataFile, MoveData, num_moves);
   	    	//printf("move data read");
   	    //	printf("MoveData akey is %lld, %lld, %lld, %lld", (long long)MoveData[0][0],(long long)MoveData[1][0],(long long)MoveData[2][0],(long long)MoveData[3][0]);
   	
/* 2.5 Read test data*/
	// here the difference is that we need to define a sampling scheme
	// proportional to the number of the new cases?
   /*double **TestData = (double**)malloc( sizeof(double *) * num_tests);
   	    	for(i = 0; i < num_tests; i++)
   	    	 {
   	    	 TestData[i] = (double*)malloc( sizeof(double) * num_tests_vars);  
   		     	}*/
   	    //@@@	read_test_data(TestDataFile, TestData, num_tests);
   	    //	printf("test data read");
   	    	
   	    	
/*2.6 Read birth data*/
    long long **BirthData = (long long**)malloc( sizeof(long long *) * num_births);
   	    	for(i = 0; i < num_births; i++)
   	    	 {
   	    	 BirthData[i] = (long long*)malloc( sizeof(long long) * num_births_vars);  
   		     	}
   	     read_birth_data(BirthDataFile, BirthData, num_births) ;
   	      // printf("birth data read");
   	      
/*2.7 Read slaughter data*/
	long long **SlaughterData = (long long**)malloc(sizeof(long long*)*num_slaughter) ;
			for(i = 0; i< num_slaughter; i++)
			{
				SlaughterData[i] = (long long*)malloc(sizeof(long long)*num_slaughter_vars) ;
			}
		   read_slaughter_data(SlaughterDataFile, SlaughterData, num_slaughter) ;
		  // printf("SlaughterData is %lld %lld %lld", SlaughterData[0][1], SlaughterData[1][1], SlaughterData[2][1]);
		  
//printf("a") ;
/*2.7 Prepare sequence data*/
if(disease==1)
{
	read_sequence_data(SequenceData,seq_length,Seq_master) ;
}
else if(disease==0)
{
	for(i=0; i < seq_length -1; i++)
	{
		Seq_master[i]= (char)"ATGC"[rand()%4];
	}
	Seq_master[seq_length-1] = '\0' ;
}

int i = 0;
 	        // while(Seq_master[i]!='\0')
        //				{
        	//				printf("%c",Seq_master[i]) ;
			//			i++ ;
			//			}

/*Create IsolateData frame that stores isolate information*/

struct isolate **head_list_isolate = (struct isolate**)malloc(sizeof(struct isolate*)*num_max_isolate) ;
//double **IsolateData = (double**)malloc( sizeof(double *) * num_max_isolate);
//char* IsolateSeqData[num_max_isolate] ; //= (char*)malloc(sizeof(char*)*num_max_isolate);
   	    	
	//printf("IsolateData created");

/*Create OutPut data*/
  //  double **OutPut = (double**)malloc( sizeof(double *) *tot_iterations);
  // 	    	for(i = 0; i < tot_iterations; i++)
   //	    	 {
   //	    	 OutPut[i] = (double*)malloc( sizeof(double) *num_OutPut); 
	//			   		     }
				   		     
/*Create FarmProductionStatus*/

      	/*============================================================================================*/
 
/* Set memory for linked lists*/      
struct animal_node **animal_node_pointer = (struct animal_node**)malloc( sizeof(struct animal_node*) * num_total_animals);	


//animal_node_pointer = (struct animal_node*)malloc( sizeof(struct animal_node) * num_total_animals);	
// animal_node_pointer points to the address of the first animal of the array with a length of num_total_animals
//struct animal_node *animal_node_pointer = (struct animal_node*)malloc( sizeof(struct animal_node) * num_total_animals);
struct animal_node *inf_slaughtered[sim_days] ; // a pointer to animal_node that stores slaughtered and infected animals
// maybe need to add slaughtered and non-infected linked list too

struct animal_node* FarmProductionList[num_farm_production]; // pointer to first animal at each farm-production
struct inf_move *Farm_ptr_inf_move[num_total_farms] = {NULL}; // this records moves that had infected animals from the source farm
struct event_node* event_day[sim_days]; // a pointer to a struct of event

// make a 3-column table that stores isolates data

int column_isolate_table = 7 ; 
//[0] id [1] rid [2] isolated time [3] indicator for E if this is (0) not chosen (1) eligible for timeframe but not chosen (2) chosen as sample
//[4] indicator for R 0 not chosen 1 chosen
//[5] indicator for H_I 0 not chosen 1 chosen
//[6] indicator for H 0 not chosen 1 chosen




int row_hid_rid = 45965 ;
int column_hid_rid = 2 ;
double **hid_rid_Data = (double**)malloc(sizeof(double*)*row_hid_rid) ;
for(i=0; i<row_hid_rid; i++)
{
	hid_rid_Data[i] = (double*)malloc(sizeof(double)*column_hid_rid) ;

}
    read_hid_rid(hid_rid, hid_rid_Data,row_hid_rid) ;
    //printf("read hid_rid") ;
    


     read_num_active_herds_rid(num_active_herds_rid,num_active_herds_rid_array,n_rid) ;
     
/*===============================================================================================================================================================================================================
------------------START OF THE SIMULATION-----------------------------------------------------------
============================================================================================*/  

int var0, var1_1,var1_2,var1_3,var2,var3,var4,var5,var6,init_sample_limit;
long random_id;

/*===========START OF ITERATION=======================================================================*/
while(iteration<tot_iterations)
//for(iteration=0; iteration<tot_iterations; iteration++)
{//iteration starts
printf("ite is %d\n",iteration);
int rid_inf_array[n_rid] = {0} ;

	for(i = 0; i<n_rid; i++)
	{
		rid_status[i]=0 ;
	//	printf("rid is %d",rid_status[i]) ;
	}
double **FarmProductionStatus = (double**)malloc( sizeof(double *) * num_farm_production);
for(i = 0; i < num_farm_production; i++)
      	  	{
     	       /*Then, for each row, specifies the required memory size for column*/
     	     FarmProductionStatus[i] = (double*)malloc( sizeof(double) * column_prostatus); 
     	      FarmProductionList[i] = NULL; // initialise the animal struct
			  /* what kind of farm variable do we need?*/ 
			  for (j=0; j < column_prostatus; j++)
			  {
			  	FarmProductionStatus[i][j] = 0 ; 
				  // initialise the status table, but keep herd id 
			  }
     	    }
     	    
     	    
double **FarmDataStatus =malloc( sizeof(double *) * num_total_farms2);
   	  	 for(i = 0; i < num_total_farms2; i++)
      	  {
     	       /*Then, for each row, specifies the required memory size for column*/
     	     FarmDataStatus[i] = malloc( sizeof(double) * num_farm_var2); 
			  /* what kind of farm variable do we need?*/ 
     	   }
//clock_t begin = clock() ;
	srand((unsigned)time(NULL));		    	
    //  printf("Iteration started");
      var0 = 0;
      died_out = &var0 ;
      random_id = rand()%num_animals;
     // int var1 = 0 ;
     // nth_isolate = &var1;
      var1_1 = 0;
      nth_isolate_detected = &var1_1 ;
      var1_2 = 0;
      nth_between = &var1_2 ; 
	  var1_3 = 0;
	  max_isolate_reached =&var1_3;
      
      var2 = 0;
      nth_slaughter = &var2;
      
      var3 = 0;
      c_transmission = &var3;
      
      var4 = 0;
      c_mutation = &var4;
      
      var5 = 0;
      c_detection = &var5;
      init_sample_limit = 0;
      sample_limit = &init_sample_limit;
      // prepare IsolateData and IsolateSeqData
     // printf("err0") ;
      for(i = 0; i < num_max_isolate; i++)
   	    	 {
   	     	head_list_isolate[i] = NULL ;
   	    	 //IsolateData[i] = (double*)malloc( sizeof(double) * num_isolate_vars);  
   	    	 //IsolateSeqData[i] = (char*)malloc(sizeof(char)*634) ;
   		    	}
      
    //  printf("err01") ;
    //  printf("nth_slaughter is %d", *nth_slaughter);
      char* OutIsolateDataFile;
      char* OutNexusFile;
	  char* OutNexusFile_E;
	  char* OutNexusFile_R;
	  char* OutNexusFile_H_I;
	  char* OutNexusFile_H;
      char number[4] ;
      // making file name for IsolateData
      snprintf(number, 4, "%d", iteration) ;
      OutIsolateDataFile = malloc(strlen(OutIsolateDataFile_origin)+4+1+4) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutIsolateDataFile, OutIsolateDataFile_origin);
      
      strcat(OutIsolateDataFile, number) ;
      strcat(OutIsolateDataFile, extension) ;
      
      OutNexusFile = malloc(strlen(OutNexusFile_origin)+4+1+4) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutNexusFile, OutNexusFile_origin);
      strcat(OutNexusFile, number) ;
      strcat(OutNexusFile, extension_nexus) ;
      
      OutNexusFile_E = malloc(strlen(OutNexusFile_origin)+4+1+5) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutNexusFile_E, OutNexusFile_origin);
      strcat(OutNexusFile_E, number) ;
      strcat(OutNexusFile_E, extension_nexus_E) ;
      
      OutNexusFile_R = malloc(strlen(OutNexusFile_origin)+4+1+5) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutNexusFile_R, OutNexusFile_origin);
      strcat(OutNexusFile_R, number) ;
      strcat(OutNexusFile_R, extension_nexus_R) ;
      
      OutNexusFile_H_I = malloc(strlen(OutNexusFile_origin)+4+1+7) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutNexusFile_H_I, OutNexusFile_origin);
      strcat(OutNexusFile_H_I, number) ;
      strcat(OutNexusFile_H_I, extension_nexus_H_I) ;
      
      OutNexusFile_H = malloc(strlen(OutNexusFile_origin)+4+1+5) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutNexusFile_H, OutNexusFile_origin);
      strcat(OutNexusFile_H, number) ;
      strcat(OutNexusFile_H, extension_nexus_H) ;
      
       for(i=0; i < num_total_animals; i++)
		{
			animal_node_pointer[i] = NULL ;
		}
      	
      	var6 = 0;
		current_id_isolate = &var6; //reset the counter for current id for isolate
      	// INITIALISE THE FARM TB DETECTED STATUS, SUM OF OCCULT, SUM OF DETECTABLE
      	for(i=0; i< num_total_farms2; i++)
      	{
      		
      		FarmDataStatus[i][0] = 0; //disease status
      		FarmDataStatus[i][1] = 0;//detected or not
      		FarmDataStatus[i][2] = 0;//R3
      	    FarmDataStatus[i][3] = 0;//R4
      	    FarmDataStatus[i][5] = 0;

      		// initialise Farm_ptr_inf_move 
      		//Farm_ptr_inf_move[i]->des_farm_id = -1 ;
      		if(disease==0)
      		{
      			random_year = 0 ;
      		if(FarmData[i][3]==2) // if STB
      		{
      			random_year = rand()%2 ; // get 0 or 1
      			
			}
			else if (FarmData[i][3]>=3) // if STD or STT
			{
				random_year = rand()%3 ;
			}
				FarmDataStatus[i][4] = random_year ;
			  }
      		
		  }
		  
				
   	    	  
   		 //   }
   		  //  printf("Isolate 0 0 is %lld", IsolateData[0][0]) ;
   		  //  system("pause") ;
   		     	
   		   //   for(i=0;i<tot_iterations;i++)
 		//	 {
  		//		for (j=0; j < num_OutPut; j++)
			//	{
			//	OutPut[i][j] = 0 ; 
			//	}
 		//	 }
      	
/* Following codes to make linked lists need to be repeated in each iteration*/
/*=====Prepare Vector of day "event_day" that links to EventNode that stores information for any events (move/born/culling, infection status change)*/

/*=====2.2.1 Farm-production pointer 1: Pointer to table of infection status=======================*/
	  	// This is table storing number of infected, susceptible and so on for each management unit
	  	

     //	printf("err02") ;
     	/*====2.2.2 Farm-production pointer 2: Pointer to animal array===========================*/   	 		  
          for(i = 0; i < sim_days; i++)
                {
                event_day[i] = NULL;
                inf_slaughtered[i] = NULL ;
                }
          // printf("err03") ;
      //	printf("event is fine");
      	/* Create a vector of size of number of animals that store a pointer to animal_node*/
      	
      	
		  //for (i=0; i < num_total_animals; i++)
      //	{
      	//	animal_node_pointer[i] = (struct animal_node*)malloc( sizeof(struct animal_node)); 
      		 //initialise
      		//Farm_ptr_inf_move[i] = NULL ;
		 // }
    //   printf("animal_node_pointer is fine"); 
        
/*======================ADDING ANIMAL NODE=============================================================*/ 		  
		  //struct animal_node *new_node;
		  struct event_node *new_event;
		  // prepare a pointer that points to animal_node and event_node 
		  long eligible_counter = 0 ;
		  long eligible_counter1 = 0;
		  long eligible_counter2 = 0;
		  long random_eligible = rand()%num_eligible+1 ;
		  long random_eligible1 = rand()%20+1 ;;
		  long random_eligible2 = rand()%20+1 ;;
		  long long initial_pro_id = -100 ;
for (i=0; i < num_animals; i++)
          { 
                
                //new_node = (struct animal_node*)malloc(sizeof( struct animal_node )); 
                long long current_akey = (long long)AnimalData[i][0];
                animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
                animal_node_pointer[current_akey]->akey = current_akey;
                animal_node_pointer[current_akey]->age_day = (int)AnimalData[i][2];
                animal_node_pointer[current_akey]->type = (int)AnimalData[i][3];
              
                
                long long current_farm = (long long)AnimalData[i][1] ; // farm id
               // int current_type = (int)AnimalData[i][3] ; // which age group (production unit)
                current_pro_id = AnimalData[i][1]*3 + AnimalData[i][3] ; // initialise the current_pro_id
                animal_node_pointer[current_akey]->current_pro_id = current_pro_id;
				animal_node_pointer[current_akey]->sex = (int)AnimalData[i][4] ;
                animal_node_pointer[current_akey]->breed = (int)AnimalData[i][5];
                animal_node_pointer[current_akey]->pregnant_status = (int)AnimalData[i][6];
			 	animal_node_pointer[current_akey]->disease_status = 0 ;
				
            	
            /*If pregnant, determine when they give birth. This data can be stored at Event table rather than the linked animal list*/
           /* IGNORE PREGNANCY FOR NOW
		    if ((int)AnimalData[i][6] ==1) // if pregnant, then add an event list to the day
            {
            	int birth_date = 70;
            	
            	
            	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
            	new_event -> akey = current_akey ;
            	
            	new_event -> production_id = current_pro_id ;
            	new_event -> event_type = 2; // 2 is giving birth
            	new_event -> next_node = NULL;
            	add_event_node(event_day,birth_date, new_event) ;
            	//visualize_list(event_day, 70);this works
            	//printf("This animal is %lld",current_akey);
          	
			}  
		
            */
        
          
	/*---------Disease status BLOCK------------------------------------------------------------------------------*/

	/*======SEED INFECTION========================================================================================*/
            // TB status - 2 ways, either extract from data or choose here
            /* Not relevant here but when updating infection pressure,
            it is easier to have an indicator if the number of occult/detectable
            are changed since previous update. Create a column called inf_chaned*/
            
                //new_node -> disease_status = (int)AnimalData[i][7];
        /* Here choose a random number below num_animals and infect this animal*/
        /* This is only choosing one infected animal.
		But it's not true because for every each animal one random number is picked so 
		there might be 0 infection but can be as many as whole by chance*/
		
        if(AnimalData[i][8]==1)
                {
                	eligible_counter++;
				}
				//oops if eligible counter reaches the random eligible here,
				//and then do not increase for a while, then it keeps seeding
				
		//infect all the animals in 39784 adults, and see if disease spreads to 41684 and 42048
		if(current_pro_id==initial_pro_id-1)
		{
			eligible_counter1++;
		}
		if(current_pro_id==initial_pro_id-2)
		{
			eligible_counter2++;
		}
			
        if ((eligible_counter == random_eligible && AnimalData[i][8]==1)||(eligible_counter1 == random_eligible1 && current_pro_id==initial_pro_id-1)|| (eligible_counter2 == random_eligible2 && current_pro_id==initial_pro_id-2))
        {
       // printf("This animal is infected %ld", random_eligible) ;
      //  printf("initial infected farm id %lld", current_pro_id) ;
          if (eligible_counter == random_eligible && AnimalData[i][8]==1)
          {
         	initial_pro_id = current_pro_id ;
	      }
      
        //	system("pause") ;
        if(disease==0)
        {
        	animal_node_pointer[current_akey]->disease_status = 2;
		}
		else if(disease==1)
		{
			animal_node_pointer[current_akey]->disease_status = 2;//change from 1 to 2 as of 26Feb
		}
        	int disease_status ;
			disease_status = animal_node_pointer[current_akey]->disease_status ;
        	//new_node -> disease_status = 1 ; // infected
        	//struct isolate *new_isolate ;
        	//new_isolate = (struct isolate*)malloc(sizeof( struct isolate )); //new_isolate is a pointer, so it's compatible to pass this to ptr_isolate
       //     printf("err05") ;
            head_list_isolate[*current_id_isolate] = (struct isolate*)malloc(sizeof(struct isolate)) ;
			head_list_isolate[*current_id_isolate]->id = *current_id_isolate ;
		//	printf("err06") ;
			//new_isolate -> id = 0 ; // this is the seed
            strcpy(head_list_isolate[*current_id_isolate]->sequence,Seq_master);
        //    printf("err07") ;
            //int a = 0;
           // while(Seq_master[a]!='\0')
        	//			{
        	//				printf("%c",Seq_master[a]) ;
			//			a++ ;
			//			}
		//				a = 0;
           //  while(new_isolate->sequence[a]!='\0')
        	//			{
        	//				printf("%c",new_isolate->sequence[a]) ;
			//			a++ ;
			//			}
            head_list_isolate[*current_id_isolate]-> parent_id = -1 ;
            head_list_isolate[*current_id_isolate]-> d_time = 0 ;
            head_list_isolate[*current_id_isolate]-> num_mut = 0 ;
            head_list_isolate[*current_id_isolate]->isolated = 0;
            head_list_isolate[*current_id_isolate]->isolated_time = -1;
            head_list_isolate[*current_id_isolate]->akey = current_akey ;
            head_list_isolate[*current_id_isolate]->pro_id = current_pro_id;//is this ok? does it change as animal moves?
            animal_node_pointer[current_akey]->ptr_isolate = head_list_isolate[*current_id_isolate];
          //  printf("err06") ;
        	//new_node -> ptr_isolate = new_isolate; //assigning a pointer new_isolate to a pointer ptr_isolate
        	FarmProductionStatus[current_pro_id][disease_status+1]++; // increment the counter for infected
        //	printf("%f animal is infected",FarmProductionStatus[current_pro_id][2] ) ;
        //	system("pause") ;
			FarmDataStatus[current_farm][0] = 1; // farm is infected
			if(disease_status>=2)
			{
				FarmDataStatus[current_farm][disease_status]++ ;
			}
        
        	(*current_id_isolate)++;
        //	printf("iso id is %d",*current_id_isolate);
        	if(*current_id_isolate>=num_max_isolate)
        	{
        		printf("max isolate reached") ;
        		system("pause") ;
			}
		} 
		else
		{
		
			//new_node -> disease_status = 0 ;
			animal_node_pointer[current_akey]->ptr_isolate = NULL ;
			//new_node -> ptr_isolate = NULL ;
			FarmProductionStatus[current_pro_id][1]++;
		}
                
            
            ///initial condition: INFECTIOUS (FMD) OR REACTIVE (TB)
            if (animal_node_pointer[current_akey]->disease_status==2)
			{
				day_to_S4 = rand()%max_S3toS4 +min_S3toS4  ; // randomly select a date to occult
				//FarmProductionStatus[current_pro_id][2] = FarmProductionStatus[current_pro_id][2]+1; // increment by 1
				
				if(day_to_S4<sim_days)
				{
				
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				
				new_event -> akey = current_akey;
			//	printf("I %lld will become occult after %d days", current_akey, day_to_occult);
				new_event -> src_pro_id = current_pro_id ;
				new_event -> des_pro_id = -100 ;
				new_event -> event_type = 7; // 6 is from S" to S3
            	new_event -> next_node = NULL;
            	new_event -> src_testarea = -100 ;
            	add_event_node(event_day,day_to_S4, new_event) ;
                   }
			}    
            
            
            //INITIAL CONDITION: EXPOSED OR OCCULT
			if (animal_node_pointer[current_akey]->disease_status==1)
			{
				//randomly select from some distributions
				day_to_S3 = rand()%max_S2toS3 +1  ; // randomly select a date to occult
				//FarmProductionStatus[current_pro_id][2] = FarmProductionStatus[current_pro_id][2]+1; // increment by 1
				
				if(day_to_S3<sim_days)
				{
				
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				
				new_event -> akey = current_akey;
			//	printf("I %lld will become occult after %d days", current_akey, day_to_occult);
				new_event -> src_pro_id = current_pro_id ;
				new_event -> des_pro_id = -100 ;
				new_event -> event_type = 6; // 6 is from S" to S3
            	new_event -> next_node = NULL;
            	new_event -> src_testarea = -100 ;
            	add_event_node(event_day,day_to_S3, new_event) ;
            //	printf("add S3");
            //	system("pause") ;
                }
			
			}    
            // finally increment for the total number of animals in this type
			FarmProductionStatus[current_pro_id][0] ++;
			
			
	/*-------------DISEASE STATUS BLOCK END--------------------------------------------------------------------*/
			
			// then initialise the next_node
			animal_node_pointer[current_akey]->next_node = NULL;
             //   new_node -> next_node = NULL; 
			 animal_node_pointer[current_akey]->previous_node = NULL ;  
                //new_node -> previous_node = NULL ;
            /* ADD THE NEW NODE TO THE ARRAY */
                 add_animal_node(FarmProductionList, current_pro_id, animal_node_pointer[current_akey] ) ;
           }  
/*======ADDING ANIMALS DONE=====================================================================================================*/      
   	
    /*
    ------------------------------------------------------------------------------
    Above section should be one-off. Reading data is only once.
    ------------------------------------------------------------------------------

// here specifies iterations

/* 3 Updating status*/
	/* Procedure.
	1) Calculate the totatl RATE of Markov events (within-herd and wildlife-to-farm.
	2) Decide a date to next Markov event based on the RATE.
	3) If this date is earlier than next Non-Markov event, then do the Markov.
	4) Else, do Non-Markov.Update the event table.
	5) Update farm status, go back to 1). */
	
	
	/*So again, do we have to track all demographic and movements on uninfected farms?
	
	At least we need to know who is providing cattle to bTB infected farms.
	
	But how about false-positive? WHen it happens, their movements will stop for a while.
	Have to model test positive even bTB status is negative.
	Then when it happens stop moving and put culling.
	*/


/*===============================================================================================
==== CREATE AND UPDATE NON-MARKOV EVENTS
================================================================================================*/
int current_event_type ;
		/*==================================ADDING MOVEMENT TO EVENT==================================================*/
		/* 1.Adding movement to the event node starts*/
		for (i=0; i < num_moves; i++)
          { 
          
                int current_day = (int)MoveData[i][1] ;
               // printf("current day is %d",current_day);
               // system("pause");
                if(current_day<sim_days)
                {
				
                current_event_type =  (int)MoveData[i][4];
                struct event_node *new_event;
                new_event = (struct event_node*)malloc(sizeof( struct event_node ));
                new_event -> event_type = current_event_type ; // 0 is calf, 1 heifer 2 adult
                long long current_akey;
				current_akey = (long long)MoveData[i][0] ;//akey
                long long current_pro_id = (long long)MoveData[i][2]*3 + current_event_type;
                
                new_event -> akey = current_akey;
                new_event -> src_pro_id = current_pro_id ;
				new_event -> des_pro_id =  (long long)MoveData[i][3]*3+current_event_type; // destination farm
                new_event -> src_testarea = (int)MoveData[i][5] ;
				new_event -> next_node = NULL;   
                /* ADD THE NEW NODE TO THE ARRAY */
                 add_event_node(event_day, current_day, new_event) ;
                 if (animal_node_pointer[current_akey]==NULL)//if this animal does not exist initially
                 {
                 	animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
                 		FarmProductionStatus[current_pro_id][0]++;
                 		FarmProductionStatus[current_pro_id][1]++; //assume that these animals are susceptible
                 		animal_node_pointer[current_akey]->akey=current_akey;
                 		animal_node_pointer[current_akey]->current_pro_id = current_pro_id ;
                 		//new_animal->breed=6;//now just assume it's unknown
                 		animal_node_pointer[current_akey]->breed = 6 ;
                 		//new_animal->sex=1; //now assume all female but if needs to be precise, have to get the data from tlogevity and add to movement data
                 		animal_node_pointer[current_akey]->sex=1;
                  		animal_node_pointer[current_akey]->type=current_event_type;
						animal_node_pointer[current_akey]->disease_status=0;
				 		animal_node_pointer[current_akey]->num_births=0;
				 		animal_node_pointer[current_akey]->type=current_event_type;
				 		animal_node_pointer[current_akey]->pregnant_status=0;
						 animal_node_pointer[current_akey]->ptr_isolate = NULL;
				 		animal_node_pointer[current_akey]->next_node=NULL;
				 		animal_node_pointer[current_akey]->previous_node=NULL;
				 		add_animal_node(FarmProductionList, current_pro_id, animal_node_pointer[current_akey] ) ;
				 }
				 
               }
           }/*Adding movement to event node is done*/
		/*==============ADD MOVEMENT DONE==============================================================*/
           
		/*==============ADD BIRTH START=================================================================*/      

		/* 2.ADD BIRTH EVENT AND ADD BIRTH ANIMALS To THE ANIMAL POINTER*/
    
          for (i=0; i < num_births; i++)
          { 
          int current_day = (int)BirthData[i][1] ;
          if(current_day<sim_days)
          {
          	long long current_akey = (long long)BirthData[i][0] ;//akey
		  if (animal_node_pointer[current_akey]==NULL)//if this animal does not exist initially
                 {
                				
                long long current_pro_id = (long long)BirthData[i][2]*3 ;
                struct event_node *new_event;
                new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
                new_event -> event_type = 3;          
                new_event -> akey = current_akey;
                new_event -> src_pro_id = -100 ;
				new_event -> des_pro_id = current_pro_id; // destination farm
                new_event -> src_testarea = -100 ;
				new_event -> next_node = NULL;          
                add_event_node(event_day, current_day, new_event) ;
                 
                //THEN CREATE ANIMAL NODE
                animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
                animal_node_pointer[current_akey]->akey = current_akey;
                animal_node_pointer[current_akey]->age_day = 0;
                animal_node_pointer[current_akey]->type = 0;
                animal_node_pointer[current_akey]->num_births = 0;
                animal_node_pointer[current_akey]->breed = (int)BirthData[i][4];
                animal_node_pointer[current_akey]->sex = (int)BirthData[i][3];
                animal_node_pointer[current_akey]->pregnant_status = 0;
                animal_node_pointer[current_akey]->disease_status=0; //important assumption that vertical transmission does not occur
                animal_node_pointer[current_akey]->current_pro_id = current_pro_id;
				animal_node_pointer[current_akey]->next_node = NULL;
                animal_node_pointer[current_akey]->previous_node = NULL;
                animal_node_pointer[current_akey]->ptr_isolate = NULL;
                //animal_node_pointer[current_akey] = new_animal; //now added memory location to the pointer list
                
            	}
            }
           }
		/*========================ADD BIRTH DONE==============================================*/

		/*=========================ADD SLAUGHTER DONE==========================================*/
 		for (i=0; i < num_slaughter; i++)
          { 
          int current_day = (int)SlaughterData[i][1] ;
          if(current_day<sim_days)
          {
          	long long current_akey = (long long)SlaughterData[i][0] ;//akey
          //	printf("akey is %lld\n",current_akey);
          //	system("pause") ;
		 // if (animal_node_pointer[current_akey]==NULL)//wait, do I need this?
          //       {
          if(current_akey <=num_total_animals)
          {
          	struct event_node *new_event;
                new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
                new_event -> event_type = 4;          
                new_event -> akey = current_akey;
                new_event -> src_pro_id = -100 ;
				new_event -> des_pro_id = -100; // destination farm
                new_event -> src_testarea = -100 ;
				new_event -> next_node = NULL;          
                add_event_node(event_day, current_day, new_event) ;
		  }
                
            //	}
            }
           }
		/*=================ADD SLAUGHTER DONE============================================================*/   
 
		/*================================ADD UNIT CHNAGE EVENTS=================================================================================*/
       int YEARS = sim_years; //@@@
           int current_day;
	    for(i=0;i<YEARS;i++)
           {
           	current_day = 365*i + 10;
           	struct event_node *new_event; // I think just defining this once is enough (i.e. no need to include this in the loop)
                new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
                
                new_event -> event_type = 8;
                new_event -> src_pro_id = -100 ; 
				new_event -> des_pro_id = -100; // destination farm
                new_event -> src_testarea = -100 ;
				new_event -> next_node = NULL;   
                 add_event_node(event_day, current_day, new_event) ;
		   }
		 /*================================ADD UNIT CHNAGE EVENTS DONE============================================================================*/
 
		/*=========================DATA PREPARATION ENDS=============================================*/ 
  

/*=========================SIMULATION===================================================*/   
double today_date = 0;
int next_non_markov_date, random_date;
double day_to_markov,updated_date;
long long move_akey, animal_akey;
long long current_akey, src_pro_id, src_farm_id,pre_move_det_status ;
double sum_inf_pressure;
struct event_node *current_event ;
struct event_node *adding_new_event ;
struct animal_node *current_animal ;

int current_pro_id,after_move_det_status,current_year,testarea_num,current_farm ;
int day_to_add ;
int j = 0;

/*==================HERE INFECT ONE OF HEIFER AND ADULT IN THE SEEDING FARM==============================*/



/*=======SIMULATION STARTS=================================================================================================================================================================*/	
while((today_date<sim_days)&& (*died_out==0) && (*max_isolate_reached==0) && (*sample_limit==0))
{
//printf("L968: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
//movesystem("pause") ;
//printf("this is day %f\n", today_date);
if(today_date>365)
					{
					
					detection_pressure_origin = detection_pressure_origin1 ;
					
					}
					else
					{
						
						detection_pressure_origin = detection_pressure_origin0 ;
					}
//add test here				
current_year = floor(today_date/365) ;
if (disease==0)
{	
if(today_date-current_year*365==0) // beginning of the year
{
	random_date = rand()%365 ;
	struct event_node *new_event; 
    new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
    new_event->akey = -1 ;
    new_event->event_type = 5 ;
    new_event->des_pro_id = -100;
    new_event->src_pro_id = -100;
    new_event->next_node = NULL;
    new_event->src_testarea = -1 ;
    //printf("adding test") ;
	add_event_node(event_day, (int)(today_date+random_date), new_event) ;
}
// add testing done	
}
	
  next_non_markov_date = ceil(today_date);

  while (event_day[next_non_markov_date] == NULL)
            {//loop2
		    next_non_markov_date++; // go to the next day
	  //  printf("go to next day") ;
            } // now get next non-markov date
  if(next_non_markov_date>=sim_days)
  {
  	break;
  }
  updated_date=update_markov_date(today_date,FarmData,FarmDataStatus,FarmProductionStatus,FarmProductionList,num_farm_production,beta_a,beta_b,next_non_markov_date, event_day,
  detection_scenario, disease, mutation_model, max_S2toS3,seq_length, current_id_isolate, mu, nth_isolate_detected,head_list_isolate, c_transmission, c_mutation, c_detection, 
  detection_pressure_origin,detection_pressure_increased, died_out,max_isolate_reached,num_isolate_upper,sample_limit,rid_status); //@need to change arguments
if(*sample_limit==1)
{
	break ; //finish simulation
}

/*==============MARKOV DID NOT HAPPEN=====================================================================*/
  if (updated_date==next_non_markov_date) // this means markov event did not happen
     {

     current_event = event_day[next_non_markov_date];
     while(current_event!=NULL)
     {
  //  printf("current event is %d", current_event->event_type) ;
     current_akey = current_event->akey ;
  //  printf("current akey is %lld\n", current_akey) ;
  //  printf("current pro id is %lld\n",current_event->src_pro_id) ;
    if (current_event-> event_type ==5 )//if this is testing
	      {
	 	//	printf("test") ;
	 	   test_farms(FarmData,FarmDataStatus,FarmProductionStatus, FarmProductionList, event_day, current_event,Se_S3, Se_S4, current_year, next_non_markov_date) ; // get the testing schedule id
	      } 
	else if(current_event->event_type==10)
		  {
		 // 	printf("lift ban") ;
		  	 src_farm_id = floor(current_event->src_pro_id/3) ;
		  	 FarmDataStatus[src_farm_id][1] = 0 ; // farm detected status goes back to 0
		  	 
		  }
     else if(animal_node_pointer[current_akey]!=NULL)
     {
     	if(animal_node_pointer[current_akey]->current_pro_id!=-10)
     	{
		 
	/*==============EVENT EITHER MOVEMENT, NEW BIRTH OR CULL/DEATH================================*/
	   if (current_event-> event_type <= 4 ) // if movement or new birth or cull death
	      {
	      //	printf("move or slaughter") ;
	        src_pro_id = animal_node_pointer[current_akey]->current_pro_id ;
	       // printf("A2") ;
	        src_farm_id = floor(src_pro_id/3) ;
	        //printf("A3") ;
	    //   printf("farm id is %lld", src_farm_id) ;
	        pre_move_det_status = FarmDataStatus[src_farm_id][1] ; //if detected 1, non detected 0
	       // printf("A4") ;
	        testarea_num = current_event->src_testarea ;
	       // printf("A5") ;
	       
			
			move_animal_unit(FarmProductionList,FarmDataStatus,FarmProductionStatus,event_day,current_event,animal_node_pointer,Se_S3,Se_S4, disease, Se_slaughter, next_non_markov_date,
		    nth_slaughter, nth_between,inf_slaughtered, head_list_isolate, Farm_ptr_inf_move, num_isolate_upper, sample_limit,FarmData, rid_status
			); // function to move animals
	
			after_move_det_status =FarmDataStatus[src_farm_id][1] ;

		  
	/*===========Here testing all animals in the herd plus contact tracing if infected animals are detected==========*/
/*	if(after_move_det_status == pre_move_det_status +1) // means this farm is detected in the movement process, either pre-movement or slaughter
	{
		// add whole herd testing schedule on this farm
		// add contact tracing
		long long inf_source_farm_id = FarmData[src_farm_id][9]; //@@@@ get who infected this farm (and when?)
	
	//can I use test_farms function to do the job?
	// how can I do test after all animals are culled? maybe add test event after one day, or 5 days?
	// contact tracing - within 5 days too?
	struct event_node *new_event; // I think just defining this once is enough (i.e. no need to include this in the loop)
    new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
    new_event -> event_type = 5;
    new_event -> src_pro_id = src_pro_id ; //this is testing schedule
	new_event -> des_pro_id = -100; // destination farm
    new_event -> src_testarea = testarea_num ;
	new_event -> next_node = NULL;   
         
     //ADD THE NEW NODE TO THE ARRAY 
    add_event_node(event_day, next_non_markov_date+5, new_event) ;
    } */
	/*==========Movement section done==========*/
	//	printf("movement done\n") ;
	      } // if movement done
	/*=====================EVENT MOVEMENT, BIRTH, CULL/DEATH DONE===========================================*/
	   
		    // if testing is happening, check all animals with test accuracy, then pick up test positive ones.
		    //Record information of positive animal (age, region of herd etc) and cull it.
		    //For positive herds, do follow up testing
	/*===============UPDATE DISEASE STATUS========================================================================*/
	   else if (current_event-> event_type == 6||current_event-> event_type == 7) //if updating disease status
	      {
	     // 	printf("event is %lld today is %d", current_event->event_type,next_non_markov_date) ;
	      	
	   //   printf("updating disease status") ;
	      animal_akey = current_event-> akey ;
	        //current_animal = animal_node_pointer[animal_akey] ;
	        //current_pro_id = current_animal->current_pro_id;
	         current_pro_id = animal_node_pointer[animal_akey]->current_pro_id ;
	         
	        current_farm = (int)(floor(current_pro_id/3)) ;
	      //  printf("current farm is %d", current_farm) ;
		    if (current_event-> event_type ==6 && animal_node_pointer[animal_akey]->disease_status == 1) //exposed to occult(TB) or exposed to infectiousness
		      {
		      	
		      //	printf("This is S2 to S3") ;
		      day_S3toS4 = rand()%(max_S3toS4)+min_S3toS4;
		 //    printf("day S3 to S4 is %d",day_S3toS4) ;
		      day_to_add = day_S3toS4 + next_non_markov_date; //define when occult to detectable or infectious to recover happens
		//	  printf("day_to_add is %d",day_to_add) ;
				  	
		      animal_node_pointer[animal_akey]-> disease_status = 2 ;
		    //  printf("Num of status 2 is %f",FarmProductionStatus[current_pro_id][2] );
		      FarmProductionStatus[current_pro_id][2]--;
		    //  printf("Num of status 3 is %f",FarmProductionStatus[current_pro_id][3] );
		      FarmProductionStatus[current_pro_id][3]++;
		      FarmDataStatus[current_farm][2]++;
		      if(FarmDataStatus[current_farm][2]!=(FarmProductionStatus[current_farm*3][3]+FarmProductionStatus[current_farm*3+1][3]+FarmProductionStatus[current_farm*3+2][3]))
		      {
		      	//printf("L1395: now FarmD is %lld and FPS is %lld - %lld - %lld",FarmData[current_farm][7],FarmProductionStatus[current_farm*3][3],FarmProductionStatus[current_farm*3+1][3],FarmProductionStatus[current_farm*3+2][3]);
		      	system("pause") ;
			  }
			  if(FarmDataStatus[current_farm][2]<=0)
              {
              //	printf("L1392:inf is %lld farm is %ld\n",FarmData[current_farm][7],current_farm) ;
              	system("pause") ;
              }
		      
		      if(day_to_add<sim_days)
			  {
		//	  	printf("will add new event for S3 to S4") ;
		      adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      adding_new_event->akey=animal_akey;
		      adding_new_event->src_pro_id = -100; 
		      adding_new_event->des_pro_id = -100;
		      adding_new_event->src_testarea=-100;
			  adding_new_event->event_type=7;//occult to detectable happens
			  
			  add_event_node(event_day, day_to_add, adding_new_event) ;
		      }
		//	  printf("S2 to S3 done") ;
		      }
		      else if (current_event-> event_type ==6 && animal_node_pointer[animal_akey]->disease_status == 3)
		      {
		 //     	printf("updated cancelled") ;
			  }
		    if (current_event-> event_type ==7 && animal_node_pointer[animal_akey]->disease_status==2) // occult to detectable (TB), infectious to immune (FMD)
		    //have to check how receovered are treated in coalescent/BD model
		      {
		      animal_node_pointer[animal_akey]-> disease_status = 3 ;
		      FarmProductionStatus[current_pro_id][3]--;
		      FarmProductionStatus[current_pro_id][4]++;
		      if(FarmProductionStatus[current_farm*3][3]==0&&FarmProductionStatus[current_farm*3][2]==0&&
		      FarmProductionStatus[current_farm*3][4]==0&&FarmProductionStatus[current_farm*3+1][4]==0&&
			  FarmProductionStatus[current_farm*3+1][3]==0&&FarmProductionStatus[current_farm*3+1][2]==0&&
			  FarmProductionStatus[current_farm*3+2][3]==0&&FarmProductionStatus[current_farm*3+2][2]==0&&
			  FarmProductionStatus[current_farm*3+2][4]==0)
		      {
		      	FarmDataStatus[current_farm][0] = 0 ;
			  }
		      FarmDataStatus[current_farm][3]++ ;//increase S4 (immune FMD or infectious TB)
		      FarmDataStatus[current_farm][2]-- ;//decrease occult
		       if(FarmDataStatus[current_farm][2]!=(FarmProductionStatus[current_farm*3][3]+FarmProductionStatus[current_farm*3+1][3]+FarmProductionStatus[current_farm*3+2][3]))
		      {
		      //	printf("L1440: now FarmD is %lld and FPS is %lld - %lld - %lld",FarmData[current_farm][7],FarmProductionStatus[current_farm*3][3],FarmProductionStatus[current_farm*3+1][3],FarmProductionStatus[current_farm*3+2][3]);
		      	system("pause") ;
			  }
		      if(FarmDataStatus[current_farm][2]<0)
		      {
		      //	printf("L1445:now it's negative %lld day %f farm is %d\n",FarmData[current_farm][7],next_non_markov_date,current_farm) ;
		      //	printf("inf is %f immune is %f pro is %lld\n",FarmProductionStatus[current_pro_id][3],FarmProductionStatus[current_pro_id][4],current_pro_id) ;
			  system("pause") ;
			  }
		      //animal_node_pointer[animal_akey]->ptr_isolate->isolated_time = next_non_markov_date;
		   			   
		      }
		      else if (current_event-> event_type ==7 && animal_node_pointer[animal_akey]->disease_status==3)
		      {
		  //    	printf("update cancelled") ;
			  }
	       // printf("Updating TB status done") ;	
	      }//exposed->occult ot occult->detectable DONE
	/*===========================UPDATE DISEASE STATUS DONE==========================================*/
	      else if (current_event->event_type==8)
	      {
	     // 	printf("Changing unit!");
	      	move_production_type(FarmProductionStatus,FarmProductionList,num_total_farms, FarmDataStatus); // @double check if I updated # of animals for each status?
	      //	printf("Changing unit done!");
	      //	system("pause");
		  }
	/*================Movement-ban lift==============================================================================*/
		 
	}

	else
	{
	//	printf("This animal is already slaughtered\n");
		//system("pause") ;
	 }
}	
		struct event_node *previous_event;
		//previous_event =  (struct event_node*)malloc(sizeof( struct event_node )); // passing a pointer that points to a memory allocated dynamically
	   previous_event = current_event ; //rewire to the next event
	   current_event = current_event->next_node;
	   if (current_event!=NULL)
	   {
	 //  	printf("going to free event now") ;
	   free(previous_event);//this is not dynamic memory, but is it ok to free?
	 //  printf("next event is %d, %lld", current_event->event_type, current_event->akey);
       }
       else
       {
     // 	printf("next event is NULL") ;
	   }
	   event_day[next_non_markov_date] = current_event;
	 
	 /*---------------------------------------------------
	 Here add capturing the number of infected and do sampling.
	 Hold on: when randomly sampling, should I include all exposed ones or only infectious ones?
	 ---------------------------------------------------------------------------------------------*/ 
	  
     }//while loop for going to next events ends
	
	
     } // this is the end of loop 1    
     today_date = updated_date;
   
  } //end of while loop
  /*=============================SIMULATION ENDS==============================================================================================================================================================*/	
	
	    
//count_farms_infected_detected(FarmData,OutPut,num_total_farms,iteration);
//printf("infected is %lf and detected is %lf",OutPut[iteration][0],OutPut[iteration][1]);
//system("pause");
//printf("wrote results");
//free memory
  /*Clear FarmPrduction*/
  //struct animal_node* node;
  //struct animal_node* temp;
   //for(i = 0; i < num_farm_production; i++)
  //       { 
  // free(FarmProductionStatus[i]); // I don't have to free this memory
  //free(FarmProductionList[i]);
  //	node = FarmProductionList[i] ;
  //	while(node!=NULL)
  //	{
  //		temp = node;
  //		free(temp) ;
  //		node = node->next_node ;
  //	}
  
   //                  }
 
   
if((*nth_isolate_detected)>=num_isolate_lower)
     {
     // make isolate table now?
	// printf("Y1") ;
	 double **isolate_table = (double**)malloc(sizeof(double*)*(*nth_isolate_detected)) ;
       for(i=0; i<*nth_isolate_detected; i++)
       {
	   isolate_table[i] = (double*)malloc(sizeof(double)*column_isolate_table) ;
       }
	   for(i=0; i<*nth_isolate_detected; i++)
          {

             for(j=0 ; j <column_isolate_table; j++)
              {
    	      isolate_table[i][j]=0;
	          }
}

    for(i=0; i<num_total_farms2; i++)
	{
	if(FarmDataStatus[i][5]==1) //if isolated
	{
		rid_inf_array[(int)FarmData[i][5]]++ ;
		}	
		
	 } 	
    //here using FarmDataStatus calculate how many farms are infected in each region
	// then get the total number of farms active as of 2002 using active data
	//(1) make an array for number of herds in each region lenhgth 16
	//(2) make an array for numebr of infected herds in each region length 16 	
     	
     	
   //  	printf("Y2") ;
    export_IsolateData(OutIsolateDataFile, OutNexusFile,OutNexusFile_E,OutNexusFile_R,OutNexusFile_H_I,OutNexusFile_H, head_list_isolate, MembershipData, num_max_isolate, nth_isolate_detected, num_isolate_upper, seq_length,current_year,isolate_table,hid_rid_Data, rid_inf_array,num_active_herds_rid_array, n_rid2) ;
	iteration++ ;
	/*Export each subsampling lists*/
	 for(j=0 ; j <*nth_isolate_detected; j++)
              {
    	      free(isolate_table[j]);
	          }
	free(isolate_table) ;
	
//	printf("Y3") ;
	 }
	 else
	 {
	 	printf("sample size is %d\n", *nth_isolate_detected) ;
	 }
   
   free(OutIsolateDataFile) ;
   
 //  printf("transmission happened %d", *c_transmission) ;
  //  printf("mutation happened %d", *c_mutation) ;
  //   printf("detection happened %d", *c_detection) ;
  //   printf("between farm tras happened %d",*nth_between) ;
   // export_seq_nexus(OutNexusFile,nth_isolate,nth_isolate_detected,seq_length,IsolateData,IsolateSeqData) ;
     free(OutNexusFile) ;
     free(OutNexusFile_E);
     free(OutNexusFile_R);
     free(OutNexusFile_H_I);
     free(OutNexusFile_H);
     
//	 printf("A") ;
	 for(i=0; i < num_total_animals; i++)
   	 {
   	 	if(animal_node_pointer[i]!=NULL)
   	 	{
   	 			animal_node_pointer[i]->ptr_isolate=NULL ; // initialise
   	 	free(animal_node_pointer[i]);
			}
   	 
		}
	//	printf("B") ; 
for(i=0; i<num_max_isolate;i++)
{
	if(head_list_isolate[i]!=NULL)
   	    	 {
   	    	 free(head_list_isolate[i]);  
   		    }
}

   		    	
   	 
	//	printf("C") ;
		struct inf_move *next_inf;
		struct inf_move *ptr;
	for(i=0; i < num_total_farms; i++)
	{
		while(Farm_ptr_inf_move[i]!=NULL)
		{
			ptr = Farm_ptr_inf_move[i] ;
			next_inf = ptr->next_node ;
			free(ptr) ;
			Farm_ptr_inf_move[i] = next_inf ;
		}
		
     	     free(FarmDataStatus[i]); 
	}
	 free(FarmDataStatus);
	
	for(i=0;i<num_farm_production;i++)
	{
		free(FarmProductionStatus[i]) ;
	}
//	clock_t end = clock() ;
//double time_spent = (double)(end-begin)/CLOCKS_PER_SEC ;
//printf("time is %f", time_spent) ;
//	printf("D") ;
free(FarmProductionStatus) ;

//and finally free event node if they are remaining
struct event_node *event;
struct event_node *next_event;
for(i=0;i<sim_days;i++)
{
	while(event_day[sim_days]!=NULL)
	{
		
		event = event_day[sim_days];
		next_event = event->next_node ;
		free(event);
		event_day[sim_days] = next_event ;
		
	}
}
}//END OF EACH ITERATION


/*===========END OF ONE ITERATION==============================================================================================================================================*/


/*=======================EXPORT OUTPUT (END OF MAIN)========================================================================================================*/
//write_OutPut(OutPutFile,OutPut,tot_iterations,num_OutPut) ;    
 for(i = 0; i < num_total_farms; i++)
      	  {
     	     free(FarmData[i]); 
     	   }
     	   free(FarmData);
     	   free(OutIsolateDataFile_origin);
     	   free(OutNexusFile_origin) ;
} //END OF MAIN INT


/*======================EXPORT OUTPUT DONE=====================================================================================================================*/
            
/*=====FUNCTIONS===============================================================================================================================*/


/* -------------------------------------------------------------------------- */
/* read_farm_data: READING AND PARSING CSV FARM LIST */
/* -------------------------------------------------------------------------- */
void read_farm_data(char FarmDataFile[], double **FarmData,int num_farms)
{     
    /* OPEN INPUT FILE */
    FILE *Farms = fopen(FarmDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    double farm_id, testarea, disease_status;
    double x_coord, y_coord, rid;
    int line_num;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_farms; line_num++)
      { 
         fscanf(Farms, "%lf,%lf,%lf,%lf,%lf, %lf",&farm_id, &x_coord, &y_coord, &testarea,&disease_status, &rid);

         /* STORE VALUES IN FARM LIST */
             FarmData[line_num][0] = farm_id;
             FarmData[line_num][1] = x_coord;
             FarmData[line_num][2] = y_coord;
             FarmData[line_num][3] = testarea ;
             FarmData[line_num][4] = disease_status; 
             FarmData[line_num][5] = rid;
             
            
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Farms);
   
} 
/*---------------------------------------------------------------------------------*/
/* READ COMMUNITY MEMBERSHIP DATA*/
/*----------------------------------------------------------------------------------*/
void read_membership_data(char MembershipDataFile[],int **MembershipData, int num_total_farms2)
{
	 /* OPEN INPUT FILE */
    FILE *Data = fopen(MembershipDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int farm_id, y0,y1,y2,y3,y4,y5,y6,y7,y8,y9 ;
    int line_num;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_total_farms2; line_num++)
      { 
         fscanf(Data, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",&farm_id, &y0, &y1, &y2,&y3,&y4,&y5,&y6,&y7,&y8,&y9);

         /* STORE VALUES IN FARM LIST */
             MembershipData[line_num][0] = y0;
             MembershipData[line_num][1] = y1;
             MembershipData[line_num][2] = y2 ;
             MembershipData[line_num][3] = y3; 
             MembershipData[line_num][4] = y4;
			 MembershipData[line_num][5] = y5; 
			 MembershipData[line_num][6] = y6;
			 MembershipData[line_num][7] = y7;
			 MembershipData[line_num][8] = y8;
			 MembershipData[line_num][9] = y9; 
			 
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Data);
	
}
/* -------------------------------------------------------------------------- */
/* READ INITIAL ANIMAL DATA*/
void read_animal_data(char AnimalDataFile[], long long **AnimalData, int num_animals)
{
	
/* OPEN INPUT FILE */
    FILE *Animals = fopen(AnimalDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, farm_id, age_day,Type, Sex, Breed, pregnant_status, disease_status, eligible;
    long long akey;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_animals; line_num++)
      { 
         fscanf(Animals, "%lld,%d,%d, %d, %d,%d,%d,%d, %d",&akey, &farm_id, &age_day, &Type, &Sex, &Breed,&pregnant_status,&disease_status, &eligible );

         /* STORE VALUES IN AnimalData */
             AnimalData[line_num][0] = akey;
             AnimalData[line_num][1] = farm_id;
             AnimalData[line_num][2] = age_day;
             AnimalData[line_num][3] = Type ;
             AnimalData[line_num][4] = Sex; 
             AnimalData[line_num][5] = Breed;
             AnimalData[line_num][6] = pregnant_status;
             AnimalData[line_num][7] = disease_status; 
             AnimalData[line_num][8] = eligible;
              
             
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Animals);
 //  printf("First animal is %lld",AnimalData[0][0]) ;
} 

/*---------------------------------------------------------------------------*/
/* READ TESTING DATA*/
void read_test_data(char TestDataFile[], double **TestData, int num_tests)
{
	
/* OPEN INPUT FILE */
    FILE *Tests = fopen(TestDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, day, test_schedule_id;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_tests; line_num++)
      { 
         fscanf(Tests, "%d,%d",&day, &test_schedule_id);

         /* STORE VALUES IN AnimalData */
             TestData[line_num][0] = day;
             TestData[line_num][1] = test_schedule_id;
             
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Tests);
   
}

 /* -------------------------------------------------------------------------- */
/* add_stub: ADD animals TO production_pointer */
/* -------------------------------------------------------------------------- */
void add_animal_node(struct animal_node *FarmProductionList[], int current_pro_id, struct animal_node *node_to_add )
{     

 struct animal_node* current_node1;
 current_node1 = FarmProductionList[current_pro_id]; // thought X[a] = *(X+a) but seems X[a] is a pointer, so it's address
 if(current_node1 == NULL)
    {
        FarmProductionList[current_pro_id] = node_to_add;
    }
 else
    {
       current_node1->previous_node = node_to_add ;
       node_to_add -> next_node = current_node1;
       FarmProductionList[current_pro_id] = node_to_add;
       
    }

}


 /* -------------------------------------------------------------------------- */
/* add_stub: ADD event TO event_day POINTER */
/* -------------------------------------------------------------------------- */
void add_event_node(struct event_node *event_day[], int day, struct event_node *node_to_add )
{     

 struct event_node *current_node1;
 // printf("event day is %p", event_day[day]);
 
 current_node1 = event_day[day];
 if(current_node1 == NULL)
    {
    	// printf("current node is null and then connect to %lld \n", node_to_add -> akey);
        event_day[day] = node_to_add;
    //    printf("this was the first event of the day %d",day);
        // printf("event day is now pointing to %p", event_day[day]);
    }
 else
    {
    //	printf("current node is %lld", current_node1 -> akey);
       node_to_add -> next_node = current_node1;
    //   	printf("next node is now %lld", node_to_add -> akey);
       event_day[day] = node_to_add;
    //    printf("current node is now %lld",event_day[day]->akey);
    //    system("pause");
    }

}


/*------------------------------------------------------------------------------
READ MOVEMENT DATA
-------------------------------------------------------------------------- */
void read_movement_data(char MoveDataFile[], long long **MoveData, int num_moves)
{     
    /* OPEN INPUT FILE */
    FILE *Moves = fopen(MoveDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    long long line_num, akey,src_farm,des_farm;
    int  day, age_type, src_testareanum;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_moves; line_num++)
      { 
         fscanf(Moves, "%lld,%d,%lld,%lld, %d, %d",&akey, &day, &src_farm, &des_farm, &age_type,&src_testareanum);

         /* STORE VALUES IN FARM LIST */
             MoveData[line_num][0] = akey;
             //printf("akey is %lld",(long long)MoveData[line_num][0]);
             //system("pause");
             //system("pause");
             MoveData[line_num][1] = day;
             MoveData[line_num][2] = src_farm;
             MoveData[line_num][3] = des_farm;
             MoveData[line_num][4] = age_type;
             MoveData[line_num][5] = src_testareanum;
      }
   /* CLOSE INPUT FILE */
   fclose(Moves);
   
} 
/*-----------------------------------------------------------------------------*/

/* READ BIRTH DATA AND MAKE EVENTS NODE*/
void read_birth_data(char BirthDataFile[], long long **BirthData, int num_births)
{     
    /* OPEN INPUT FILE */
    FILE *Births = fopen(BirthDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    long long line_num, akey;
    int bdate, farm, sex, breed;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_births; line_num++)
      { 
         fscanf(Births, "%lld,%d,%d,%d,%d",&akey, &bdate, &farm, &sex, &breed);
         
             BirthData[line_num][0] = akey;
             BirthData[line_num][1] = bdate;
             BirthData[line_num][2] = farm;
             BirthData[line_num][3] = sex;
             BirthData[line_num][4] = breed;
      }
      fclose(Births);
}

/*--------------READ SLAUGHTER DATA----------------------------------------*/
void read_slaughter_data(char SlaughterDataFile[], long long **SlaughterData, int num_slaughter)
{
	/* OPEN INPUT FILE */
    FILE *Slaughter = fopen(SlaughterDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    long long line_num, akey;
    int cull_date ;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_slaughter; line_num++)
      { 
         fscanf(Slaughter, "%lld,%d",&akey, &cull_date);
         
             SlaughterData[line_num][0] = akey;
             SlaughterData[line_num][1] = cull_date;
             
      }
      fclose(Slaughter);
	
	
}
/*--------------------------------------------------------------------------------*/
/*-------------READ hid_rid DATA------------------------------------------------------*/
void read_hid_rid(char hid_rid[], double **hid_rid_Data, int row_hid_rid)
{
	FILE *Temp = fopen(hid_rid,"r") ;
	int line_num;
	double farm_id, rid;
	for(line_num = 0; line_num < row_hid_rid ; line_num++)
	{
		fscanf(Temp, "%lf, %lf", &farm_id, &rid) ;
		hid_rid_Data[line_num][0] = farm_id ;
		hid_rid_Data[line_num][1] = rid ;
	}
	fclose(Temp) ;
}

/*-----------READ num_active_herds_rid-------------------------------*/
void read_num_active_herds_rid(char num_active_herds_rid[], int num_active_herds_rid_array[], int row_active)
{
	FILE *Temp = fopen(num_active_herds_rid,"r") ;
	int line_num,n ;
	for(line_num = 0; line_num < row_active ; line_num++)
	{
		fscanf(Temp, "%d", &n) ;
		num_active_herds_rid_array[line_num] = n ;
		
	}
	fclose(Temp) ;
}



/*---------------------------------------------------------------------------------
-- Random seed
----------------------------------------------------------------------------------*/
/*
int random_farm_seed(double P, int num_farms,double **FarmData, double **FarmProductionStatus)
{
	int num_p = roundf(P*num_farms) ;
	
	int i=0, j;
	int array[num_p];
	int current_farm;
	int current_herd_size;
	while(i<num_p)
	{
		array[i] = rand()%(num_farms-1);
		for(j=0;j<i;j++)
		if(array[j]==array[i])
		break;
		if(!(j<i))
		{
		current_farm = array[i];
		FarmData[current_farm][8] = 1;
		int rand_num = rand()%99;
		int indicate;
		if(rand_num<33)
		{
			indicate = 0;//infected is calf unit
		}
		else if (rand < 67)
		{
			indicate = 1;//infected is heifer unit
		}
		else 
		{
		indicate = 2;//infected id adult unit - they may not have the equal probability of infection
		}
		FarmProductionStatus[current_farm+indicate][1] = 1; //infect this management unit
		// so now allocated infection to FarmData and FarmProductionStatus
		// then we have to select infected animal on the infected premises
		// this will be perhaps easier to do on the linked list
		i++;
		}
		
	}
	
}
*/
/*=============================================================================
-- 
*/ 
/* -------------------------------------------------------------------------- */
/* VISUALIZE INFORMATION FROM LINKED LIST ----------------------------------- */
/* -------------------------------------------------------------------------- */
void visualize_list(struct event_node *event_day[], int day)  
{

 
 struct event_node *current_node1;
 current_node1 = event_day[day];
 //printf("current node is %lld",current_node1->akey);
 if(current_node1 != NULL)
    {
    //   printf("Day %d: ", day );
       while(current_node1 != NULL)
          {
          //	if (current_node1->des_pro_id==(1225*3+2)||current_node1->src_pro_id==(1225*3+2))
           //   {
		//	 printf("%d,%lld ", current_node1 -> event_type, current_node1 -> akey);
			 // system("pause") ;
		   // }
              
              current_node1 = current_node1 -> next_node; 
              
          }  
      // printf("\n");
   }
   
 

}

/* -------------------------------------------------------------------------- */
/* VISUALIZE LINKED ANIMALS ----------------------------------- */
/* -------------------------------------------------------------------------- */
void visualize_animals(struct animal_node *FarmProductionList[], int production_id)  
{

 int nth = 0;
 struct animal_node *current_node1;
 current_node1 = FarmProductionList[production_id];
 if(current_node1 != NULL)
    {
     //  printf("Id is %d: ", production_id);
       while(current_node1 != NULL)
          {
          //   if(current_node1->akey==2454078)
           //  {
			 
		//	 printf("nth %d current animal is %lld \n",nth, current_node1 -> akey);
             
             // system("pause") ;
         // }
              current_node1 = current_node1 -> next_node; 
              nth ++ ;
              
              
          }  
       //printf("\n");
   }
   
 

}

/* -------------------------------------------------------------------------- */
/* MOVE ANIMAL NODE FROM ONE TO OTHER FARM */
/* -------------------------------------------------------------------------- */
void move_animal_unit(struct animal_node *FarmProductionList[],double **FarmDataStatus,double **FarmProductionStatus,struct event_node *event_day[], struct event_node *current_event, 
struct animal_node **animal_node_pointer,double Se_occult, double Se_detect, int disease, double Se_slaughter, int date, int *nth_slaughter, int *nth_between,
struct animal_node *inf_slaughtered[], struct isolate **head_list_isolate, struct inf_move *Farm_ptr_inf_move[], int num_isolate_upper, int *sample_limit, double **FarmData, int rid_status[])
{
	int today_date = date;
	int current_event_type = current_event->event_type;
	//printf("event is %d", current_event_type) ;
	long long current_akey = current_event->akey ;
	double random_value;
	long long des_farm_id,src_farm_id ;
	int stop = 0;
	struct animal_node *moving_animal;
	long long src_pro_id ;
	int disease_status ;
if(animal_node_pointer[current_akey]==NULL)
	{
		printf("this animal is null and event is %d",current_event_type);
		current_event_type = -10 ; // ignore this event
		system("pause") ;
	}
else //else animal_node_pointer[] is not NULL
{
	src_pro_id = animal_node_pointer[current_akey]->current_pro_id ;
	src_farm_id = (long)floor(src_pro_id/3) ;
	moving_animal = animal_node_pointer[current_akey] ;
	disease_status = moving_animal->disease_status ;

}	//else animal_node_pointer[] is not NULL DONE


	
 /*========Movement to slaughter==============================================================================================================*/	
 if(current_event_type==4) //slaughter/death
 {
 //	printf("slaughter event");
 	//system("pause") ;
   
 	if(src_pro_id!=-10)
 	{
 /*	if(src_pro_id==69196)
 	{
 		printf("L1696: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
        //system("pause") ;
	 }*/
 	//change FarmProductionStatus values
 //	printf("pro id is %lld",src_pro_id) ;
 	FarmProductionStatus[src_pro_id][0]--;
 	FarmProductionStatus[src_pro_id][disease_status+1]--; //minus 1
 	if(disease_status>=2) //if TB and the animal is REACTIVE OR INFECTIOUS: No slaughter inspection for FMD.
 	{
 	//	printf("Farm det is %lld", FarmData[src_farm_id][6]) ;
 		if(disease==0 && FarmDataStatus[src_farm_id][1]==0)//if TB and not detected yet
 		{
 		//	printf("almost detection") ;
 		random_value = (double)rand()/(double)RAND_MAX ;
 			if(random_value<=Se_slaughter)
 			{
 			FarmDataStatus[src_farm_id][1]=1; // becomes detected 
 			if(FarmDataStatus[src_farm_id][5]==0)
 			{
 			moving_animal->ptr_isolate->isolated = 1 ;
			moving_animal->ptr_isolate->isolated_time =(double)today_date;
			(*nth_isolate_detected)++ ;
			FarmDataStatus[src_farm_id][5]=1;	
		             	if(rid_status[(int)FarmData[src_farm_id][5]]==0)
						   {
						   rid_status[(int)FarmData[src_farm_id][5]]=1 ;
						   }
			 }
 			
			if(*nth_isolate_detected>=num_isolate_upper)
			{
				*sample_limit = 1 ;
				printf("simulation stop") ;
				return;
			}
			
			//printf("nth is %d time is %f",*nth_isolate_detected,moving_animal->ptr_isolate->isolated_time) ;
			// add movement ban lift
				if((today_date+movement_ban_days)<sim_days)
			  	{
			  		int day = today_date+movement_ban_days ;
			  	        struct event_node *adding_new_event;
			  			adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      			adding_new_event->akey=-1; // akey does not matter
		      			adding_new_event->src_pro_id = src_pro_id; 
		      			adding_new_event->des_pro_id = -100;
		      			adding_new_event->src_testarea=-100;
			  			adding_new_event->event_type=10;//movement ban lift
			         //   printf("set up adding node") ;
			  
			            struct event_node *current_node1;
                        current_node1 = event_day[day];
                        if(current_node1 == NULL)
                         {
                         event_day[day] = adding_new_event;
						 }
                         else
                         {
                         adding_new_event -> next_node = current_node1;
                         event_day[day] = adding_new_event;
                         }
                      //   printf("node added") ;
			 			 
				  }
		 	}	
		 }
		
 		//FarmData[src_farm_id][4]=1 ; //becomes infected (should be marked as infected already)
 		FarmDataStatus[src_farm_id][disease_status] --; //minus 1
 		if(FarmDataStatus[src_farm_id][disease_status]<0)
			  			{
			  			//	printf("2133:now it's negative %lld farm %lld",FarmData[src_farm_id][disease_status+5],src_farm_id) ;
						  system("pause") ;
						  }
		if((FarmDataStatus[src_farm_id][2]+FarmDataStatus[src_farm_id][3]+FarmProductionStatus[src_farm_id*3][2]+FarmProductionStatus[src_farm_id*3+1][2]+FarmProductionStatus[src_farm_id*3+2][2])==0)
		{
			FarmDataStatus[src_farm_id][0]=0 ;
		}
						  
		
	}
/*================Remove animals from the linked lists=============================================================*/
if (moving_animal -> previous_node != NULL)
	{ 
//	printf("A");
	// if moving animal's previous node is conencted to other animal
	    //    
		struct animal_node *prev_animal;
        prev_animal = moving_animal -> previous_node ; // get the node which moving one is connected
	    if (moving_animal -> next_node != NULL) // if moving one's next node is conencted
	            {
	              //	printf("B starts");
	  	        struct animal_node *next_animal;
                next_animal = moving_animal -> next_node ; // get the next animal
                prev_animal -> next_node = next_animal;    // reconenct previous_animal's next node
                next_animal -> previous_node = prev_animal;//similarly reconnect next_animal's previous node
		         //  printf("B ends");
				}	
		else // if next node is null
		        {
		   	    prev_animal -> next_node = NULL ;
		        }
		    //  printf("A ends");    
	}
else // if previous node is null 
	{
	      //   printf("C starts");	
	    if (moving_animal -> next_node != NULL) //and if next node is not null
	    {  
	  //  printf("D starts") ;
	  	struct animal_node *next_animal;
        next_animal = moving_animal -> next_node ;
     //   printf("next animal akey is %lld",next_animal->akey) ;
        next_animal -> previous_node = NULL;
        FarmProductionList[src_pro_id] = next_animal;
      //  printf("D ends") ;
		}
		else 
		{
	//		printf("E starts") ;
		FarmProductionList[src_pro_id] = NULL ; 
		}
	//	       printf("C ends");  
	}
 moving_animal->current_pro_id = -10 ; //means it's sent to slaughterhouse 
 moving_animal->previous_node = NULL;
 moving_animal->next_node = NULL ;
/*===============Removing the animal from the linked list done=================================================================================*/		
 
	 	(*nth_slaughter)++ ;
	/*if(src_pro_id==69196)
 	{
 		printf("L1768: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
       // system("pause") ;
	 } */
 }
 //printf("slaughter ends") ;
}
/*===========================Movement to slaughter END====================================================================================*/

/*===========================Movement except new calves born===============================================================================*/	
  else if (current_event_type != 3&& (animal_node_pointer[current_akey]->current_pro_id!=-10)&& (animal_node_pointer[current_akey]->disease_status!=3)) // if not 3 it means it is not new born
  {
  	/*
  	if(src_pro_id==69196)
 	{
 		printf("L1780: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
        //system("pause") ;
	 }*/
   // printf("Move starts\n");  
//	int src_farm_id = floor(src_pro_id/3) ;
//	printf("current_event_type is %d", current_event_type);
/* Check if movement is allowed*/	
	if ((current_event_type<=2 && FarmDataStatus[src_farm_id][1]==1)||(src_pro_id!=current_event->src_pro_id)) 
	// if Disease status is detected here, cancel
	// or if src_pro_id stored in the event node is not matching to the pro_id stored in animal_node_pointer, cancel the move
	{
		stop = 1;
	//	printf("stop farm is detected");
	//	system("pause");
	}
	else if ((current_event_type==1||current_event_type==2) && current_event->src_testarea == 0 && disease ==0) //pre-movement testing if disease is TB, farm in MCA, non-calf movement
	{
	/*==========Pre-movement test=====================================================================*/
		
		int num_occult = FarmProductionStatus[src_pro_id][3] ;
		int num_detectable = FarmProductionStatus[src_pro_id][4] ;
		if (num_occult+num_detectable>0)
		{
			double P_miss = (pow(1-Se_occult,num_occult))*(pow(1-Se_detect,num_detectable)) ;
			if ( ((double)rand()/(double)RAND_MAX) > P_miss)
			/// needs to complete
			{
				FarmDataStatus[src_farm_id][1]=1 ; //farm becomes detected
				// now omit the function which animals to be detected
				//in future to add the detected animals to "detected list"
				stop = 1 ;
			//	printf("stop! pre movement detected");
			//	system("pause");
			}
			
		}
	}
	//printf("movement stop is %d",stop);
	/*===========Pre-movement test done===================================================================*/
/*Check if movement is allowed ENDS*/
	else
	{
		stop = 0;
	}


/*===================IF MOVEMENTS ALLOWED=================================================================================*/
if (stop == 0) // if the movement is still allowed
{
//printf("non-stop\n") ;
des_farm_id = floor((current_event->des_pro_id)/3);
//printf("pro id %lld", current_event->des_pro_id) ;
  	
  	/*===========Add this move to Farm_ptr_inf_move if the source farm is infected and if this move is not yet recorded*/
  	if(FarmDataStatus[src_farm_id][0]==1) // if this farm is infected
  	{
  		if((disease==1&&(disease_status==1||disease_status==2))||(disease==0&&(disease_status>=1)))
  		{
  	
  			(*nth_between)++;
  		//	printf("infected moves %d",*nth_between) ;
  			
  			
		  }
  	if(disease_status==2||disease_status==3)
  	{
  		FarmDataStatus[src_farm_id][disease_status]--;
  			if(FarmDataStatus[src_farm_id][disease_status]<0)
  			{
  				printf("L2343: negative!") ;
  				system("pause") ;
			  }
			 
	  }
  		//system("pause") ;
  	//	printf("This %lld farm is infected",src_farm_id) ;
/*
  		if(FarmData[des_farm_id][9]==-1)
  		{
  			FarmData[des_farm_id][9] = src_farm_id; // record who infected this farm
		  }
  		if(Farm_ptr_inf_move[src_farm_id]==NULL)
  		{
  		//	printf("A\n") ;
  		//	printf("src_farm_id is NULL") ;
  		struct inf_move *new_node = (struct inf_move*)malloc(sizeof(struct inf_move));
		  new_node->date = today_date ;
		  new_node->des_farm_id = des_farm_id ;
		  new_node->next_node = NULL;
		  Farm_ptr_inf_move[src_farm_id] = new_node;
		  }
		else 
		{
		//	printf("B\n") ;
	//========NEEDS CHANGE IF INFECTED TO BE RECORDED============================//		
		//	printf("src_farm_id %lld is not NULL", src_farm_id) ;
			
  // printf("The value of des is : %lld\n", Farm_ptr_inf_move[src_farm_id]->des_farm_id  );
  // system("pause") ;
		struct inf_move *next_node ;
		next_node = Farm_ptr_inf_move[src_farm_id];
			if((next_node->date!=today_date)&&(next_node->des_farm_id!=des_farm_id))
			// this if: only adding movement if destination and date is different. This needs to be changed if I want to record this movement contained infected animals.
			// So that I can go through all animals to check their disease status and record new_node.infected = 1, which is currently not one of elements in this node.
			{
				struct inf_move *new_node = (struct inf_move*)malloc(sizeof(struct inf_move));
		 	 	new_node->date = today_date ;
		 	 	new_node->des_farm_id = des_farm_id ;
		  		new_node->next_node = NULL;
				Farm_ptr_inf_move[src_farm_id] = new_node ;
				new_node->next_node = next_node;
			//	printf("new node added which is %lld", new_node->des_farm_id) ;
			//	system("pause") ;
			}
		} */
	//========NEEDS CHANGE IF INFECTED TO BE RECORDED DONE=======================
	  } // recording infected movements done
  	       
	if (moving_animal -> previous_node != NULL)
	          { // if moving animal's previous node is conencted to other animal
	    //      printf("A starts");
	  //  printf("C\n") ;
	 //  printf("previous node is %lld",moving_animal -> previous_node->akey ) ;
			  struct animal_node *prev_animal;
              prev_animal = moving_animal -> previous_node ; // get the node which moving one is connected
	            if (moving_animal -> next_node != NULL) // if moving one's next node is conencted
	               {
	    //           	printf("B starts");
	  	           struct animal_node *next_animal;
                   next_animal = moving_animal -> next_node ; // get the next animal
                   prev_animal -> next_node = next_animal;    // reconenct previous_animal's next node
                   next_animal -> previous_node = prev_animal;//similarly reconnect next_animal's previous node
		//           printf("B ends");
				   }	
		        else // if next node is null
		           {
		   	       prev_animal -> next_node = NULL ;
		           }
		//       printf("A ends");    
	           }
	else // if previous node is null 
	           {
	    //       	printf("D\n") ;
	        //  	printf("previous node is NULL") ;
	    //       printf("C starts");	
	           if (moving_animal -> next_node != NULL) //and if next node is not null
	              {  
	  	          struct animal_node *next_animal;
                  next_animal = moving_animal -> next_node ;
                  next_animal -> previous_node = NULL;
                  FarmProductionList[src_pro_id] = next_animal;
		          }
		       else 
		          {
		   	      FarmProductionList[src_pro_id] = NULL ; 
		          }
		//        printf("C ends");  
		        }	
       //}//LOOPA
       
	if (FarmProductionStatus[src_pro_id][0]>0)
	   {
	 //  	printf("E\n") ;
	   FarmProductionStatus[src_pro_id][0] = FarmProductionStatus[src_pro_id][0] - 1;
	   //int current_status = moving_animal->disease_status ;
	   FarmProductionStatus[src_pro_id][disease_status+1] = FarmProductionStatus[src_pro_id][disease_status+1] - 1;	
	    if(FarmDataStatus[src_farm_id][2]!=(FarmProductionStatus[src_farm_id*3][3]+FarmProductionStatus[src_farm_id*3+1][3]+FarmProductionStatus[src_farm_id*3+2][3]))
		      {
		      	//printf("L2484: now FarmD is %lld and FPS is %f - %f - %f\n",FarmData[src_farm_id][7],FarmProductionStatus[src_farm_id*3][3],FarmProductionStatus[src_farm_id*3+1][3],FarmProductionStatus[src_farm_id*3+2][3]);
		      //	printf("L2485: status is %d pro id is %lld farm status is %f\n", disease_status, src_pro_id,FarmData[src_farm_id][4]) ;
				  system("pause") ;
			  }
	if(disease==1) //FMD
	{
	if(FarmDataStatus[src_farm_id][0]==1&&FarmProductionStatus[src_farm_id*3][3]==0&&FarmProductionStatus[src_farm_id*3][2]==0&&
			  FarmProductionStatus[src_farm_id*3+1][3]==0&&FarmProductionStatus[src_farm_id*3+1][2]==0&&
			  FarmProductionStatus[src_farm_id*3+2][3]==0&&FarmProductionStatus[src_farm_id*3+2][2]==0)
		      {
		      	FarmDataStatus[src_farm_id][0] = 0 ;
		     // 	printf("infection clear! %lld\n", src_farm_id) ;
			  }	
	}
	else if(disease==0) //TB
	if(FarmDataStatus[src_farm_id][0]==1&&FarmProductionStatus[src_farm_id*3][3]==0&&FarmProductionStatus[src_farm_id*3][2]==0&&
			  FarmProductionStatus[src_farm_id*3][4]==0&&FarmProductionStatus[src_farm_id*3+1][3]==0&&
			  FarmProductionStatus[src_farm_id*3+1][2]==0&&FarmProductionStatus[src_farm_id*3+1][4]==0&&
			  FarmProductionStatus[src_farm_id*3+2][3]==0&&FarmProductionStatus[src_farm_id*3+2][2]==0&&
			  FarmProductionStatus[src_farm_id*3+2][4]==0)
		      {
		      	FarmDataStatus[src_farm_id][0] = 0 ;
		     // 	printf("infection clear! %lld\n", src_farm_id) ;
			  }
	    
	   
	 //and if there are no more TB animal change farm status to 0
	   
	   }
    
    

	
//Part B: add moving animals to its new farm
// this applies to no matter what the animal is

// how should I consider about new farm?
int des_pro_id = current_event->des_pro_id;
    if (FarmProductionList[des_pro_id] == NULL)
    {
    	
    	FarmProductionList[des_pro_id] = moving_animal;
    	moving_animal -> previous_node = NULL;
    	moving_animal -> next_node = NULL;
	}
	else
	{
		
	   struct animal_node *new_next_animal ;
	   
	   	new_next_animal = FarmProductionList[des_pro_id] ;
	   	new_next_animal -> previous_node = moving_animal;
	   	moving_animal -> next_node = new_next_animal;
	   	moving_animal -> previous_node = NULL ;
	   	FarmProductionList[des_pro_id] = moving_animal ;
	   	
	}
	FarmProductionStatus[des_pro_id][0]++;
	//int current_status = moving_animal->disease_status ;
	FarmProductionStatus[des_pro_id][disease_status+1]++;

	if(disease_status>=2)
	{
		FarmDataStatus[des_farm_id][disease_status]++;

	}
	if(FarmDataStatus[des_farm_id][0]==0)
	{
		if(disease==0&&disease_status>=1)
		{
			FarmDataStatus[des_farm_id][0] = 1 ;
		//	printf("now farm infected") ;
		}
		else if(disease==1 && (disease_status==1||disease_status==2))
		{
			FarmDataStatus[des_farm_id][0] = 1 ;
		//	printf("now farm infected") ;
		}
	}
		

	if(FarmDataStatus[des_farm_id][2]!=(FarmProductionStatus[des_farm_id*3][3]+FarmProductionStatus[des_farm_id*3+1][3]+FarmProductionStatus[des_farm_id*3+2][3]))
		      {
		      	
		      	//printf("L2539: now FarmD is %lld and FPS is %f - %f - %f",FarmData[des_farm_id][7],FarmProductionStatus[des_farm_id*3][3],FarmProductionStatus[des_farm_id*3+1][3],FarmProductionStatus[des_farm_id*3+2][3]);
		      	system("pause") ;
			  }
			 
//	if(current_status >= 1 && current_status<=2)
//	{
//		printf("transmission to new pro id %lld and inf is now %f",des_pro_id,FarmProductionStatus[des_pro_id][current_status+1]) ;
//	system("pause") ;
//	}
	
//	if(moving_animal==NULL)
//	{
//		printf("why null?") ;
//	}
	moving_animal->current_pro_id=des_pro_id;
} // if stop==0 ENDS

//printf("movement done\n") ;
/*==========================IF MOVEMENTS ALLOWED ENDS===========================================================*/
/*if(src_pro_id==69196)
 	{
 		printf("L2028: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
        //system("pause") ;
	 }*/
} //if (current_event_type != 3) ENDS

else if (current_event_type == 3 &&(animal_node_pointer[current_akey]->current_pro_id!=-10)) // if this is new born just add
{
	/*if(src_pro_id==69196)
 	{
 		printf("L2037: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
        //system("pause") ;
	 }*/
//	printf("event is %d",current_event_type) ;
int des_pro_id = current_event->des_pro_id;
    if (FarmProductionList[des_pro_id] == NULL)
    {
    	FarmProductionList[des_pro_id] = moving_animal;
    	moving_animal -> previous_node = NULL;
    	moving_animal -> next_node = NULL;
	}
	else
	{
	   struct animal_node *new_next_animal ;
	   
	   	new_next_animal = FarmProductionList[des_pro_id] ;
	   	new_next_animal -> previous_node = moving_animal;
	   	moving_animal -> next_node = new_next_animal;
	   	moving_animal -> previous_node = NULL ;
	   	FarmProductionList[des_pro_id] = moving_animal ;
	}
	FarmProductionStatus[des_pro_id][0]++;
	int current_status = moving_animal->disease_status;
	FarmProductionStatus[des_pro_id][current_status+1]++;
	if(current_status>=2)
	{
		FarmDataStatus[(int)floor(des_pro_id/3)][current_status]++;
		if(FarmDataStatus[(int)floor(des_pro_id/3)][0]==0)
		{
			FarmDataStatus[(int)floor(des_pro_id/3)][0]= 1;
		}
		}	
	moving_animal->current_pro_id=des_pro_id;
} // new born ENDS
//printf("J\n") ;
/*
if(src_pro_id==69196)
 	{
 		printf("L2074: T %f S %f E %f I %f", FarmProductionStatus[69196][0],FarmProductionStatus[69196][1],FarmProductionStatus[69196][2],FarmProductionStatus[69196][3]) ;
        //system("pause") ;
	 }*/
//	 printf("move done") ;
}
    
    
/* -------------------------------------------------------------------------- */
/* TESTING HERDS */
/* -------------------------------------------------------------------------- */
void test_farms(double **FarmData,double **FarmDataStatus, double **FarmProductionStatus, struct animal_node *FarmProductionList[], struct event_node *event_day[], struct event_node *current_event,double Se_S3, double Se_S4, int current_year, double today_date)
{
	int random_year, num_S3,num_S4 ;
	struct animal_node *current_animal ;
//	printf("test starts") ;
	for (i=0; i<num_total_farms; i++)
	{
		if(FarmDataStatus[i][4]==current_year) // if the farm has the same testing id
		{
			int current_testarea = FarmData[i][3] ;
			if(current_testarea>=3)
			{
				num_S3 = FarmProductionStatus[i*3+2][3] ;
				num_S4 = FarmProductionStatus[i*3+2][4] ;
				random_year = rand()%3+1 ;
			}
			else if(current_testarea==2)
			{
				num_S3 = FarmProductionStatus[i*3+2][3] ;
				num_S4 = FarmProductionStatus[i*3+2][4] ;
				random_year = rand()%2+1 ;
			}
			else if(current_testarea==1||current_testarea==0)
			{
				num_S3 = FarmProductionStatus[i*3+2][3] + FarmProductionStatus[i*3+1][3] ;
				num_S4 = FarmProductionStatus[i*3+2][4] + FarmProductionStatus[i*3+1][4];
				random_year = 1 ;
			}
			else 
			{
				printf("This farm does not have an eligible testarea!!") ;
			}
			
		  
		  if((num_S3+num_S4)>0) // if this farm has more than 0 occult+detectable
		  {
		     double P_miss = (pow(1-Se_S3,num_S3))*(pow(1-Se_S4,num_S4)) ;
			 if ( ((double)rand()/(double)RAND_MAX) > P_miss)
			  {
			//  	printf("detection") ;
			  //FarmData[i][6]=1 ; //wait until it gets slaughtered
				// now have to change the disease status of every animal
				if(FarmProductionStatus[i*3+2][3]>0||FarmProductionStatus[i*3+2][4]>0)
				{
					current_animal = FarmProductionList[i*3+2] ;
					while(current_animal!=NULL)
					{
						if(current_animal->disease_status>=2)
						{
					// change day_to_add to today_date + 1
					//include today_date in function		
						if(today_date+1<sim_days)
			  			{
			  			int day = today_date + 1 ;
			  			struct event_node *adding_new_event ;
		      			adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      			adding_new_event->akey=current_animal->akey;
		      			adding_new_event->src_pro_id = -100; 
		      			adding_new_event->des_pro_id = -100;
		      			adding_new_event->src_testarea=-100;
			  			adding_new_event->event_type=4;//slaughter
			  
			  			 struct event_node *current_node1;
                        current_node1 = event_day[day];
                        if(current_node1 == NULL)
                         {
                         event_day[day] = adding_new_event;
						 }
                         else
                         {
                         adding_new_event -> next_node = current_node1;
                         event_day[day] = adding_new_event;
                         }
		      			}
						}
						current_animal = current_animal->next_node ;
					}
				}
				if((FarmProductionStatus[i*3+1][3]>0||FarmProductionStatus[i*3+1][4]>0)&&current_testarea<=1)
				{
				current_animal = FarmProductionList[i*3+1] ;
					while(current_animal!=NULL)
					{
						if(current_animal->disease_status>=2)
						{
					// change day_to_add to today_date + 1
					//include today_date in function		
						if(today_date+1<sim_days)
			  			{
			  			struct event_node *adding_new_event ;
		      			adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      			adding_new_event->akey=current_animal->akey;
		      			adding_new_event->src_pro_id = -100; 
		      			adding_new_event->des_pro_id = -100;
		      			adding_new_event->src_testarea=-100;
			  			adding_new_event->event_type=4;//slaughter
			  
			            int day = today_date + 1 ;
			  			 struct event_node *current_node1;
                        current_node1 = event_day[day];
                        if(current_node1 == NULL)
                         {
                         event_day[day] = adding_new_event;
						 }
                         else
                         {
                         adding_new_event -> next_node = current_node1;
                         event_day[day] = adding_new_event;
                         }
		      			}
						}
						current_animal = current_animal->next_node ;
					}
				}
			  }
		  }
		   FarmDataStatus[i][4] = FarmDataStatus[i][4] + random_year ; // updating testing year
		}
	}
//	printf("test done") ;
}
/*----------------------------------------------------------------------------------*/
/* MOVING ANIMALS BETEWEEN PRODUCTION TYPE*/
/*----------------------------------------------------------------------------------*/
void move_production_type(double **FarmProductionStatus,struct animal_node *FarmProductionList[],int num_total_farms2,double **FarmDataStatus )
{
struct animal_node *moving_animal;
struct animal_node *next_animal;
struct animal_node *target_animal;
struct animal_node *previous_animal;
struct animal_node *ptr;


int current_status;
//printf("move unit");
 for (i=0; i<num_total_farms;i++)
     {
    	//printf("%d",i);
    	ptr = FarmProductionList[i*3+1] ;
    	/*if(i==6432)
    	{
    		printf("L2128: T %f S %f E %f I %f",FarmProductionStatus[i*3+1][0],FarmProductionStatus[i*3+1][1],FarmProductionStatus[i*3+1][2],FarmProductionStatus[i*3+1][3]) ;
    		printf("L2129: T %f S %f E %f I %f",FarmProductionStatus[i*3+2][0],FarmProductionStatus[i*3+2][1],FarmProductionStatus[i*3+2][2],FarmProductionStatus[i*3+2][3]) ;
		}*/
     	while(ptr!=NULL)//if there are heifer
     	{
		 if(ptr->type==1)
		 {
		 	moving_animal = ptr ;
		 	(moving_animal->type)++; 
			(moving_animal->current_pro_id)++;
		 //adding the next animal to farm pointer
		 
		
		 	if(moving_animal->next_node!=NULL) //means moving_animal is not the tail
		 	{
		 	next_animal = moving_animal->next_node ;
		 		if(moving_animal->previous_node==NULL) // means moving_animal is the top head
		 			{
		 			next_animal->previous_node=NULL;
		 			FarmProductionList[i*3+1] = next_animal;
			 		}
			 	else // means there are animals between moving_animal and FPL
			 		{
			 		previous_animal = moving_animal->previous_node;
			 		previous_animal->next_node = next_animal ;
			 		next_animal->previous_node = previous_animal ;
			 		}
		 	
	     	} //end: moving_animal is not the tail
	     	else //moving_animal is the tail
	     	{
	     		if(moving_animal->previous_node==NULL) //means moving_animal is head
		 			{
		 				FarmProductionList[i*3+1] = NULL;
		 			}
	     		else //means moving_animal is not the head
	     		{
	     			previous_animal = moving_animal->previous_node;
	     			previous_animal->next_node = NULL ;
				 }
		 	}
		 	moving_animal->previous_node = NULL;
		 	moving_animal->next_node = NULL;
		 
		 //update number of animals in the unit
		 	current_status = moving_animal->disease_status;
		 	FarmProductionStatus[i*3+1][current_status+1]--;
		 	FarmProductionStatus[i*3+1][0]--;
		 
		 
		 
		 //adding the current animal to the destination
		  // first animal in adult list
		 	if(FarmProductionList[i*3+2]==NULL)
		 	{
		 	FarmProductionList[i*3+2] = moving_animal;
		 	}
		 	else
		 	{
		 	target_animal = FarmProductionList[i*3+2];
		 	target_animal->previous_node = moving_animal;
		 	FarmProductionList[i*3+2] = moving_animal;
		 	moving_animal->next_node = target_animal;
		 	}
		 
     	 	FarmProductionStatus[i*3+2][current_status+1]++;
		 	FarmProductionStatus[i*3+2][0]++;
		 	if(FarmDataStatus[i][2]!=(FarmProductionStatus[i*3][3]+FarmProductionStatus[i*3+1][3]+FarmProductionStatus[i*3+2][3]))
		      {
		      //	printf("L2746: now FarmD is %lld and FPS is %f - %f - %f",FarmData[i][7],FarmProductionStatus[i*3][3],FarmProductionStatus[i*3+1][3],FarmProductionStatus[i*3+2][3]);
		      	system("pause") ;
			  }
		}
		 	ptr = ptr->next_node ; //move onto the next
		
	    }//while ends
	    /*	if(i==6432)
    	{
    		printf("L2206: T %f S %f E %f I %f",FarmProductionStatus[i*3+1][0],FarmProductionStatus[i*3+1][1],FarmProductionStatus[i*3+1][2],FarmProductionStatus[i*3+1][3]) ;
    		printf("L2207: T %f S %f E %f I %f",FarmProductionStatus[i*3+2][0],FarmProductionStatus[i*3+2][1],FarmProductionStatus[i*3+2][2],FarmProductionStatus[i*3+2][3]) ;
		} */
  //   	printf("heifer ends");
     	//next calf
     	ptr = FarmProductionList[i*3] ;
     	 /*	if(i==6432)
    	{
    		printf("L2214: T %f S %f E %f I %f",FarmProductionStatus[i*3+1][0],FarmProductionStatus[i*3+1][1],FarmProductionStatus[i*3+1][2],FarmProductionStatus[i*3+1][3]) ;
    		printf("L2215: T %f S %f E %f I %f",FarmProductionStatus[i*3][0],FarmProductionStatus[i*3][1],FarmProductionStatus[i*3][2],FarmProductionStatus[i*3][3]) ;
		}*/
     	while(ptr!=NULL)//if there are heifer
     	{
		 if(ptr->type==0)
		 {
		 	moving_animal = ptr ;
		 	(moving_animal->type)++;
			(moving_animal->current_pro_id)++;
		

		 //adding the next animal to farm pointer
		 
		 if(moving_animal->next_node!=NULL) //means moving_animal is not the tail
		 
		 {
		 	next_animal = moving_animal->next_node ;
		 		if(moving_animal->previous_node==NULL) // means moving_animal is the top head
		 			{
		 			next_animal->previous_node=NULL;
		 			FarmProductionList[i*3] = next_animal;
			 		}
			 	else // means there are animals between moving_animal and FPL
			 		{
			 		previous_animal = moving_animal->previous_node;
			 		previous_animal->next_node = next_animal ;
			 		next_animal->previous_node = previous_animal ;
			 		}
		 	
	     	} //end: moving_animal is not the tail	
		 	else //moving_animal is the tail
	     	{
	     		if(moving_animal->previous_node==NULL) //means moving_animal is head
		 			{
		 				FarmProductionList[i*3] = NULL;
		 			}
	     		else //means moving_animal is not the head
	     		{
	     			previous_animal = moving_animal->previous_node;
	     			previous_animal->next_node = NULL ;
				 }
		 	}
		 	moving_animal->previous_node = NULL;
		 	moving_animal->next_node = NULL;
		 
		 //update number of animals in the unit
		 current_status = moving_animal->disease_status;
		 FarmProductionStatus[i*3][current_status+1]--;
		 FarmProductionStatus[i*3][0]--;
		 
		
		 
		 //adding the current animal to the destination
		  // first animal in adult list
		 if(FarmProductionList[i*3+1]==NULL)
		 {
		 	FarmProductionList[i*3+1] = moving_animal;
		 }
		 else
		 {
		 	target_animal = FarmProductionList[i*3+1];
		 	target_animal->previous_node = moving_animal;
		 	FarmProductionList[i*3+1] = moving_animal;
		 	moving_animal->next_node = target_animal;
		 }
		 
     	 FarmProductionStatus[i*3+1][current_status+1]++;
		 FarmProductionStatus[i*3+1][0]++;
		 if(FarmDataStatus[i][2]!=(FarmProductionStatus[i*3][3]+FarmProductionStatus[i*3+1][3]+FarmProductionStatus[i*3+2][3]))
		      {
		      //	printf("L2860: now FarmD is %lld and FPS is %f - %f - %f",FarmData[i][7],FarmProductionStatus[i*3][3],FarmProductionStatus[i*3+1][3],FarmProductionStatus[i*3+2][3]);
		      	system("pause") ;
			  }
		 }
		 ptr = ptr->next_node ;
		 //if checking age, do here in future
		
		 
	    }//calf ends
	     /*	if(i==6432)
    	{
    		printf("L2295: T %f S %f E %f I %f",FarmProductionStatus[i*3+1][0],FarmProductionStatus[i*3+1][1],FarmProductionStatus[i*3+1][2],FarmProductionStatus[i*3+1][3]) ;
    		printf("L2296: T %f S %f E %f I %f",FarmProductionStatus[i*3][0],FarmProductionStatus[i*3][1],FarmProductionStatus[i*3][2],FarmProductionStatus[i*3][3]) ;
		}*/

	 }  	
//	printf("move unit done") ;	
}




/*----------------------------------------------------------------------------------*/
/* CALCULATE NEXT MARKOV EVENT DATE*/
/*----------------------------------------------------------------------------------*/
double update_markov_date(double today_date, double **FarmData,double **FarmDataStatus, double **FarmProductionStatus,struct animal_node *FarmProductionList[],int num_farm_production,double beta_a,double beta_b, 
int next_non_markov_date, struct event_node *event_day[], int detection_scenario, int disease, int mutation_model, int max_S2toS3, int seq_length, 
int *current_id_isolate, double mu, int * nth_isolate_detected,struct isolate **head_list_isolate,  int *c_transmission, int *c_mutation, int *c_detection, double detection_pressure_origin,double detection_pressure_increased, int *died_out,
int *max_isolate_reached, int num_isolate_upper, int *sample_limit,int rid_status[])
{

double inf_pressure, mutation_pressure,inf_pressure_wild,sum_pressure,cumu_pressure, detection_pressure;
sum_pressure = 0;
cumu_pressure = 0;
detection_pressure = 0;
double day_to_markov,random_value2;
int farm_id, rid ;
int k = 0;
int testareanum ;
int inf_counter = 0 ;
double tb_detection_p;
tb_detection_p=0.9; //@@@@do I put higher sensitivity for this?

//printf("update markov") ;

 for(i=0; i<num_farm_production; i++)
	{//for loop A
	
	farm_id = (int)floor(i/3) ;
	testareanum = FarmData[farm_id][3] ;
	
	/*=================Calculate the total rate of event in Gillespie========================================================*/
	if ((testareanum==0||FarmProductionStatus[i][3]>0||FarmProductionStatus[i][4]>0||FarmProductionStatus[i][2]>0)&&FarmProductionStatus[i][0]>0)
	{    //either: (farm is in movement control || there is occult animal || there is detectable animal) AND there is susceptible animal AND there is animal in the herd
	
	    if(testareanum==0&&disease==0) //if it is in MCA AND TB
	    {
	    	
	      inf_pressure_wild = beta_b*FarmProductionStatus[i][1] ; // by setting beta_b = 0, wildlife infection becomes nonexistent
	    //  printf("inf pressure wild is %f",inf_pressure_wild);
	      
		}
		else
		{
		inf_pressure_wild = 0;	
		}
		// disease specific 
		if(disease==0) //if TB 
		{
		inf_pressure = beta_a*FarmProductionStatus[i][1]*(FarmProductionStatus[i][4])/FarmProductionStatus[i][0] + inf_pressure_wild ; // bSI/N + wildlife, only [4] detectable (infectious) assumed to be infectious, not reactive stage
		if(inf_pressure<0)
		{
			printf("Error:inf_pressure cannot be negative!") ;
			system("pause") ;
		}
			
			if(mutation_model==1) 
			{
				mutation_pressure = mu*(FarmProductionStatus[i][2]+FarmProductionStatus[i][3]+FarmProductionStatus[i][4]) ; //mutation occurs in exposed, occult, and detectable equally
				FarmProductionStatus[i][7] = mutation_pressure;
				if(mutation_pressure<0)
				{
					printf("Error:mutation_pressure cannot be negative!") ;
					system("pause") ;
					}
			} //@ set mu and mutation_pressure parameter
	    }
	    else if (disease==1)//FMD
	    {
	    inf_pressure = beta_a*FarmProductionStatus[i][1]*(FarmProductionStatus[i][3])/FarmProductionStatus[i][0]  ; //bSI/N	
	  //  if((FarmProductionStatus[i][3]+FarmProductionStatus[i][2])>0)
	  //  {
	   // inf_counter ++ ;
	  //  printf("inf_counter is %d\n", inf_counter) ;
	   // if(inf_counter>1)
	   // {
	   // system("pause") ;	
	//	}
	    if(inf_pressure<0)
		{
			printf("Error:inf_pressure cannot be negative!") ;
			system("pause") ;
		}
	
		//}
		
	    if(mutation_model==1) //@ add mutation_model parameter
			{
				mutation_pressure = mu*(FarmProductionStatus[i][2]+FarmProductionStatus[i][3]); //mutation occurs in exposed and infectious equally
		//		if(mutation_pressure>0)
		//		{
		//		printf("mut pressure is %lf",mutation_pressure) ;
	   	//		system("pause") ;
		//		}
				
				FarmProductionStatus[i][7] = mutation_pressure;
			}
		}
		FarmProductionStatus[i][5] = inf_pressure; //inf_pressure, stored in column 5, indicates the total infection pressure in the herd
	/*========Calculate detection pressure and add it to the total rate=============================================*/
		if(mutation_model==1) // this needs to be changed in future if detection pressure is considered when it's not mutation model
		{
			detection_pressure = 0 ;
			// add here detection/control options
			if(FarmProductionStatus[i][3]>0||FarmProductionStatus[i][4]>0) // if TB occult/detectable, if FMD infectious and recovered
			
			{
	/* @@@@@ This part needs to be completed later: but now slaughter house inspection is implemented and this section for TB seems unnecessary			
			if(disease==0) //TB
			{
				
			if(detection_scenario==0) // every infected has the same probability
			    {
			detection_pressure = 0.24/365; // this is the daily rate of detection at slaughterhouse. Crude estimate based on 43 out of 178 reactors were found to have bTB lesions in 2015.
			// I am using this parameter for an infected herd to be detected at slaughterhouse, but the original figure is the proportion of individual animals that were reactor and found visible legions in slaughterhouse.	
			    }
			else if(detection_scenario==1@@@) // infecteds have different probabilities to be infected (or different periods to be detected: how can I consider this? switch on and off based on year?)
			// or each herd has different test timing? 
			    {
			    	detection_pressure = @@@;
				//depending on testareanum and year add or not detection_pressure in each herd
				}
			else if(detection_scenario==2@@@) // contact tracing
			    {
			    	detection_pressure = @@@;
			    	// once farm is detected, farms that were associated with will be traced back and forward.
			    	// now have to think do I know in this simulation where they came from?
				}
		    } //TB done, next FMD
	------------------------------------------------------------------------------------------------------------------@@@@*/
		    
		    
		    if(disease==1) // if FMD and FMD infectious exists yet detected - but don't care about recovery as you can't sequence from recoverd animals 
		    {
		    	if(detection_scenario==0) // every infected has the same probability
			    {
			    if(FarmDataStatus[farm_id][1]==0&&FarmProductionStatus[i][3]>0) // this farm should not be already detected
				{
					
							detection_pressure = detection_pressure_origin;
				}
			
		        }
				else if(detection_scenario==1)
				{
				  if(FarmDataStatus[farm_id][1]==0&&FarmProductionStatus[i][3]>0) //detection only if it is infectious(=clinical)
				  {
						
					//	printf("check1");
				   rid = (int)FarmData[farm_id][5] ;
				 //  printf("rid is %d",rid);	
				 //  printf("0 is %d",rid_status[0]) ;
				 //  printf("status is %d",rid_status[rid]) ;
			        if(rid_status[rid]==1)
				    {
					detection_pressure = detection_pressure_increased ;
					}	
				    else if(rid_status[rid]==0)
				    {
					detection_pressure = detection_pressure_origin ;
				    }
			      }
			// printf("check done") ;
				}
			
											 // probability of infections to be detected for FMD?! Should not be too high to allow the spread. use Bolivian parameter 0.3
											// I am using this parameter for an infected herd to be detected at passive surveillance by vets (I believe) not by farmer	
		    }
		/*@@@@ Needs to complete other scenarios below: different surveillance pressure between regions and contact tracing.	    
			else if(detection_scenario==1@@@) // infecteds have different probabilities to be infected (or different periods to be detected: how can I consider this? switch on and off based on year?)
			// what's the justification of this for FMD? maybe high risk area people are more vigilant?
			    {
			    	detection_pressure = @@@;
				//depending on testareanum and year add or not detection_pressure in each herd
				}
			else if(detection_scenario==2@@@) // contact tracing
			    {
			    	detection_pressure = @@@; 
			    	// once farm is detected, farms that were associated with will be traced back and forward.
			    	// now have to think do I know in this simulation where they came from?
			    	// how to do this for FMD? maybe just stick to forward tracing, anyway because the purpose is to see how biased samples would affect the inferneces
				} @@@@ --------------------------------------------*/
			}
				
			 // not store detection pressure but can easily calculated by [6] - [7] - [5]: CORRECTION - THIS IS NOT TRUE. [6] is a cumulative pressure over herds up to the current herd, so this formula won't give detection_pressure.
			// Should I store detection_pressure or calculate it again given that it sounds useless to store such a simple number in the data?
			// for now calculate it again but can be changed later
			
/*==============================Calculate detection pressure done=========================================================*/		
			cumu_pressure = cumu_pressure + inf_pressure + mutation_pressure + detection_pressure; // cumu_pressure is the total rate up to this herd, so this will be used when selecting random herd. So this is recorded in column 6.
			sum_pressure = sum_pressure + inf_pressure+ mutation_pressure + detection_pressure;	// sum_pressure is the total rate for ALL the herd. So this won't be recorded one by one. 
			
		}
		else // for now ignore non-mutation model: but detection_pressure should be considered in non-mutation models too when the effect of control strategy is to be studied.
		{
			cumu_pressure = cumu_pressure + inf_pressure; // calculate the cumulative pressure
			sum_pressure = sum_pressure + inf_pressure ;
		}
		//printf("sum is %lf inf is %lf mut is %lf det is %lf",cumu_pressure,inf_pressure,mutation_pressure,detection_pressure) ;
		 
		FarmProductionStatus[i][6] = cumu_pressure;
		
		
		
		//printf("inf pressure is %f",sum_pressure);
	}
	else
		{
			FarmProductionStatus[i][6] = 0 ;
		}

	/*=================Calculating total event rate done=======================================================================*/
	
	}//END OF for loop A
//	printf("end of the loop") ;
//	if(inf_counter>0)
//	{
//	printf("%d infected farm exists\n",inf_counter) ;
	//system("pause") ;	
//		}	
	
	if(sum_pressure==0)
	{
		*died_out = 1 ;
		printf("yesy dies out") ;
	//	system("pause") ;
		today_date = next_non_markov_date;
		return(today_date) ;
	}
	//printf("inf pressure is %f",sum_pressure);
	 double random_value = (double)(rand()+1)/((double)(RAND_MAX)+1);
//	 printf("random value is %f",random_value);
     day_to_markov =(-log(random_value))/sum_pressure ; // Waiting time
     if(day_to_markov<0)
     {
     	printf("Error: markov date cannot be negative!") ;
     	system("pause") ;
	 }
//  printf("day_to_markov is %f",day_to_markov);
//  printf("random is %f", random_value) ;
    // system("pause") ;
//ASSESS IF NEXT MARKOV TIME IS BEFORE NEXT NON MARKOV TIME     
	if (next_non_markov_date>day_to_markov+today_date)
	{ // if markov comes first, choose markov events 
	   k++;
	   //printf("This is %d th markov",k);
	   today_date = day_to_markov+today_date ;
	   //printf("ceil days is %d",today_date);
	   //system("pause");
	   double random_value = 0;
	   while(random_value==0)
	   {
	   	random_value =  (double)rand()/(double)(RAND_MAX)*sum_pressure;
	   }
	  
	   //if this is 0 that's trouble, which can happen
/*=======Define the farm on which an event will occur====================*/
	   
	   
	   double current_value = 0;
	   long pro_id = 0;
	   
	   while(random_value>current_value) //random_value defines which farm event will occur at
	   {
	   	
	   	if(mutation_model==1)
	   		{
	   		if(FarmProductionStatus[pro_id][6]>0) //either infection or mutation pressure has to be >0 for the farm to have events
	   		
			   {
			   current_value = FarmProductionStatus[pro_id][6] ;
	//   		printf("current value is %f",current_value) ;
	//   		printf("random_value is %f",random_value);
	//   		printf("sum_pressure is %f",sum_pressure);
	   	        }
		    }
		else if(mutation_model==0)
			{
			if(FarmProductionStatus[pro_id][5]>0)
			current_value = FarmProductionStatus[pro_id][6] ;
			}
		if(current_value>=random_value) // means the farm that has sum_pressure closest to but not bigger than random_value will be selected as the farm where an event occur
	       {
//       	printf("farm determined") ;
			   break;
	    		   }
		 else
		 {
		 	pro_id++;	
	//	 printf("pro id is %d", pro_id);
		 }    
		   //(disease==0&&(FarmData[farm_id][3]==0||FarmProductionStatus[pro_id][3]>0||FarmProductionStatus[pro_id][4]>0)&&(FarmProductionStatus[pro_id][0]!=0)) //the last one, it used to be FarmProductionStatus[pro_id][0] means if this is equal to 0, which should not happen
		   //FarmData[farm_id][3] is testarea, so checking if this is MCA
		   //next one evaluates if occult exists, next one for detectable exists
		   //ok, now think about FMD
	     
	     //  	printf("current_value is %f",current_value);
	      // 	printf("random value is %f",random_value);
	       
	   
	   }//pro_id is the farm to choose for the event
	   
	   farm_id = (int)floor(pro_id/3) ; // where farm_id was used??
	//   printf("farm id is %d", farm_id) ;
/*=======END: Define the farm on which an event will occur====================*/


/*==========================================================================================================
DEFINE WHICH EVENTS OCCURRED IN WHICH ANIMALS
(1) CHOOSE MUTATION, INFECTION OR DETECTION BASED ON RATE
(2) IF INFECTION, ALL SUSCEPTIBLE HAVE EQUAL RISKS = SIMPLY 1/S
(3.1) IF MUTATION, ALL INFECTED HAVE EQUAL RISKS 
(3.2) IF MUTATION, ALL NUCLEOTIDE HAVE EQUAL RISK
(3.3) IF MUTATION, ALL BASES HAVE EQUAL RISK TO BE TRANSFORMED TO
(4.1) IF DETECTION,
(4.2) FOR FMD ALL INFECTED ANIMALS ARE ASSUMED TO BE DETECTED, WHICH SAMPLES SHOULD I TAKE? MAYBE RANDOMLY ONE ASSUMING THAT MUTATION IS LITTLE WITHIN HERD?
(4.3) FOR BTB, TRY BOTH CASES WHERE ALL INFECTED CAN BE DETECTED OR SOME INFECTED REMAIN In THE HERD DUE TO LACK OF SENSITIVITY?
===========================================================================================================*/
	if(mutation_model==1)
	{
	//printf("mut model") ;
	// here to calculate the detection pressure again
		if(disease==1)
		{
		if(detection_scenario==0)
		    {
		    if(FarmDataStatus[farm_id][1]==0&&FarmProductionStatus[pro_id][3]>0)
			{
			detection_pressure = detection_pressure_origin ;
			}
			else
			{
			detection_pressure = 0;
			}
	//	printf("detection_pressure is %f",detection_pressure) ;
		    }
		else if(detection_scenario==1)
		    {
		    if(FarmDataStatus[farm_id][1]==0&&FarmProductionStatus[pro_id][3]>0)
			{
			   if(rid_status[(int)FarmData[farm_id][5]]==1)
			   {
			   detection_pressure = detection_pressure_increased ;
		       }
		       else
		       {
		       	detection_pressure = detection_pressure_origin ;
			   }
			}
			else
			{
			detection_pressure = 0;
			}	
		    }
	    }
	/* @@@@ detection pressure for other scenarios to be defined.-----------------------------	
		else if(detection_scenario==1)
		{
		detection_pressure = @@@;
		}
		else if(detection_scenario==2)
		{
		detection_pressure = @@@ ;
		}
	@@@@------------------------------------------------------*/	
	random_value2 = ((double)(rand()+1)/((double)(RAND_MAX)+1))*(FarmProductionStatus[pro_id][7]+FarmProductionStatus[pro_id][5]+detection_pressure); // get random value between 0 and sum of inf/mutation/detection pressure
 //  printf("random value2 is %f\n",random_value2) ;
   if(random_value2==0)
   {
   	//printf("inf is %lld mut is %lld",FarmProductionStatus[pro_id][5],FarmProductionStatus[pro_id][7]) ;
   }
	}
    else
    {
    random_value2 = (double)(rand()+1)/((double)(RAND_MAX)+1)*(FarmProductionStatus[pro_id][5]);
	}
/*=====================CHOOSE WHICH EVENTS?==============================================================*/	
//printf("random value 2 is %f", random_value) ;
/*===================================WITHIN HERD TRANSMISSION PART STARTS===========================================================================================*/
	if(random_value2<=FarmProductionStatus[pro_id][5])
	{
//	printf("This is transmission random value is %f\n",random_value2) ;
	//	printf("inf pressure in this farm is %f\n",FarmProductionStatus[pro_id][5]) ;
//	printf("c trans is %d", *c_transmission);
		(*c_transmission)++ ;
//	printf("c trans is %d\n", *c_transmission);
//	printf("c_transmission is %d",*c_transmission) ;
	//	system("pause") ;
	//	system("pause") ;
	//now get the first animal in the production unit
	   struct animal_node *current_animal;
//	 printf("yes1") ;
	   //current_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
	   if(FarmProductionList[pro_id]==NULL)
	   {
	   	printf("FPS is NULL why?") ;
	   	system("pause") ;
	   }
	   current_animal = FarmProductionList[pro_id];
//	   printf("pro id is %ld", pro_id) ;
	   int total_s = FarmProductionStatus[pro_id][1]; // number of susceptible	
	   int random_n = rand()%total_s+1;
 //printf("total_s is %d",total_s) ;
	   int counter, counter2;
	   counter = 0;
	   
	   while(current_animal!=NULL)
	    {//choosing which animal to infect begins
//	    printf("let's choose infected one'") ;
	    if(current_animal->disease_status==0) // if susceptible
	    {
	    counter++;	
//	    printf("counter is %d", counter) ;
		}
		if(counter==random_n)
			{ 
//			printf("random_n met") ;
			// add passing sequence here
			current_animal->disease_status=1;
			FarmProductionStatus[pro_id][1]--; //susceptible - 1
			FarmProductionStatus[pro_id][2]++; //Exposed + 1
	//		printf("new infected animal is determined") ;
			break;
	    	}
	    else
	     {
	     current_animal= current_animal->next_node ;		
		 }
	    
	    } 
		if(current_animal==NULL)
		{
		//	printf("animal can't be found. total is %d, current %d\n",total_s,counter) ;
		//	printf("pro id is %lld",pro_id) ;
			system("pause") ;
		}
		
		// picked up an individual that will be infected
	//    printf("while loop done") ;
	/*=======Need to decide who infected this one so that seq can be passed to===========*/
	// get random number between 0 and (number of infected animals in this herd - 1)
		int num_infectious;
		if(disease==0)
		{
		num_infectious = FarmProductionStatus[pro_id][4] ; // for TB this is the number of detectable(infectious)
		}
		else if(disease==1)
		{
		num_infectious = FarmProductionStatus[pro_id][3] ; // for FMD this is the number of infectious
		}
		int random_n2 = rand()%num_infectious + 1; // a pseudo-random number between 1 and num_infectious
		counter2 = 0 ;
		struct animal_node *infecting_animal ;
		infecting_animal = FarmProductionList[pro_id];
	
	while(infecting_animal!=NULL)
	{
		if((infecting_animal->disease_status==2&&disease==1)||(infecting_animal->disease_status==3&&disease==0))// FMD and this animal is infectious
		{
			counter2++;
		}
		if(counter2==random_n2)
		{
			head_list_isolate[*current_id_isolate] = (struct isolate*)malloc(sizeof(struct isolate)) ;
			strcpy(head_list_isolate[*current_id_isolate]->sequence, infecting_animal->ptr_isolate->sequence) ;// copy sequence
		
			head_list_isolate[*current_id_isolate]->d_time = today_date;
			head_list_isolate[*current_id_isolate]->isolated_time = -1;
			head_list_isolate[*current_id_isolate]->num_mut = 0;
			head_list_isolate[*current_id_isolate]->pro_id = pro_id;
			head_list_isolate[*current_id_isolate]->isolated = 0;
			head_list_isolate[*current_id_isolate]->parent_id = infecting_animal->ptr_isolate->id ;
			head_list_isolate[*current_id_isolate]->akey = current_animal->akey ;
		//	printf("B") ;
			 // increment isolate id before assigning: this should be done for other occasions to be consistent
			head_list_isolate[*current_id_isolate]->id = *current_id_isolate; 
			//printf("err4")	;
			
		//	printf("C") ;
			current_animal->ptr_isolate = head_list_isolate[*current_id_isolate] ; // assign the pointer
		//	printf("err5")	;
		//printf("New infection is %lld",infecting_animal->akey) ;
		//	system("pause") ;
		(*current_id_isolate)++;
		if(*current_id_isolate>=num_max_isolate)
        	{
        		*max_isolate_reached = 1 ;
			}
			break;	
		}
		else
			{
			infecting_animal= infecting_animal->next_node ;	
			}
	}
	if(infecting_animal==NULL)
	{
		printf("Infecting animal can't be found! num_inf is %d, current %d",num_infectious,counter2);
		system("pause") ;
	}
			//need to set it's event to be occult
	   	 	struct event_node *adding_new_event;
	   	 	
	   	 	int day_to_S2toS3 = rand()%max_S2toS3+1;
		    int day_to_add = day_to_S2toS3 + floor(today_date);
			  
		    if(day_to_add<sim_days)
			  {
		      adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      adding_new_event->akey=current_animal->akey;
		      adding_new_event->src_pro_id = -100; 
		      adding_new_event->des_pro_id = -100;
		      adding_new_event->src_testarea=-100;
			  adding_new_event->event_type=6;//occult to detectable happens
			  
			   struct event_node *current_node1;
                        current_node1 = event_day[day_to_add];
                        if(current_node1 == NULL)
                         {
                         event_day[day_to_add] = adding_new_event;
						 }
                         else
                         {
                         adding_new_event -> next_node = current_node1;
                         event_day[day_to_add] = adding_new_event;
                         }
		      }
	// printf("transmission done\n") ;
       	}
	      
	   
	/*===================================TRANSMISSION PART ENDS==============================================================================================*/


	else if((random_value2>FarmProductionStatus[pro_id][5])&&(random_value2<=(FarmProductionStatus[pro_id][5]+FarmProductionStatus[pro_id][7]))) // if mutation is selected
	/*==============MUTATION PART STARTS=================================================================================================================*/
		{
//	printf("This is mutation") ;
		//	system("pause") ;
	//	printf("c mut is %d", *c_mutation) ;
		(*c_mutation)++ ;
	//	printf("c mut is %d\n", *c_mutation) ;
			//get random value n which is euql to the number of suscptible or that are subjcted to mutations
			//get nth animals that are eligible
		int num_tot_mutatable = (int)(FarmProductionStatus[pro_id][7]/mu) ; //total number of animals that are subjected to mutations
		int random_n3 = rand()%num_tot_mutatable+1; //get random number between 1 and num_tot_mutatable
	//	printf("A") ;
		 struct animal_node *current_animal;
	//	 printf("B") ;
		 current_animal = FarmProductionList[pro_id];
	//	 printf("C") ;
		 if(current_animal==NULL)
		 {
		 	printf("Something wrong! Animals should exist here but not!");
		 	system("pause") ;
		 }
		 else
		 {
		 //	printf("current animal is not null") ;
		 	int counter3;
		 	counter3 = 0;
		 	while(counter3<random_n3)
		 	{
		 	if((disease==0&&current_animal->disease_status!=0)||(disease==1&&(current_animal->disease_status==1||current_animal->disease_status==2))) //TB disease==0, FMD 1. if TB then if animals are either exposed/
			 counter3++; // increment
			 
			 if(counter3==random_n3)
			 {
			 	struct isolate *mutating_isolate ;
			 	mutating_isolate = current_animal->ptr_isolate; //now have to choose which base to change
			 	/*============================================================================================
			 	MUTATION MODEL HERE*/
			 	int r_position = rand()%(seq_length-1) ; //random position at nucleotide seq (random number generator needs improvement)
				char next_base; // make sure it is just char, not char[], which becomes character string type
				char current_base ; //get the current base
    			current_base = mutating_isolate->sequence[r_position] ;
    	//		printf("current base is %c", current_base) ;
    			if(current_base=='A') // make sure the character is bracketted with '', not " "
						{
		//					printf("Matched!");
							next_base= (char)"TGC"[rand()%3];
	
						}
				else if(current_base=='T')
						{
		//				printf("Matched!");
						next_base= (char)"AGC"[rand()%3];	
						}
				else if(current_base=='G')
						{
		//				printf("Matched!");
						next_base= (char)"ATC"[rand()%3];
						}
				else if(current_base=='C')
						{
		//				printf("Matched!");
						next_base= (char)"ATG"[rand()%3];	
						}
				else
						{
						printf("Error! This is not eligible base!") ;
						system("pause") ;
						}
		//		printf("Next base is %c", next_base);
				mutating_isolate->sequence[r_position] = next_base;
		//		printf("Now base is %c", mutating_isolate->sequence[r_position]);
				(mutating_isolate->num_mut)++; //mutation count+
				//@ have to add times of mutation here
				/* MUTATION MODEL ENDS
				======================================================================================================*/
				break;
			}//ENDS of IF
			else
			{
			current_animal = current_animal->next_node ; //go to next animal	
			}	
			} // END OF WHILE
		 } //else done
	//	 printf("mut done") ;
		} // END OF MUTATION
	/*===================================MUTATION PART ENDS================================================================================================*/
	
	/*===================================DETECTION PART STARTS=======================================================================*/
	
	else //else detection
	{
	//	printf("This is detection") ;
	//	printf("detection pressure is %f", detection_pressure) ;
	//	printf("threshold is %f and ran val is %f",FarmProductionStatus[pro_id][5]+FarmProductionStatus[pro_id][7], random_value2) ;
	//	system("pause") ;
//	printf("c det is %d", *c_detection) ;
		(*c_detection)++ ;
	//	printf("c det is %d", *c_detection) ;
		struct animal_node *current_animal;
		 current_animal = FarmProductionList[pro_id];
		 struct animal_node *previous_animal;
		 struct animal_node *next_animal ;
		 struct animal_node *next_slaughtered ;
		 
		 //struct isolate *current_isolate; 
		 if(current_animal==NULL)
		 {
		 	printf("Something wrong! Animals should exist here but not!");
		 	system("pause") ;
		 }
		 if(FarmDataStatus[farm_id][1]==1)
		 {
		 	printf("This farm is already detected and should not be re-detected again!");
		 	system("pause") ;
		 }
		 else //else1
		 {
		 //	printf("yep") ;
		 //	printf("This farm has animals and not detected yet") ;
		 	long num_inf = FarmDataStatus[farm_id][2] ;//number of infectious
		 	if(num_inf<0)
		 	{
		 		printf("Error: num inf cannot be negative!") ;
		 		system("pause") ;
			 }
		 //	printf("num inf is %lld", num_inf) ;
	//	 printf("num inf is %ld", num_inf) ;
		 	int random_inf = rand()%num_inf; // get random number between 0 and num_inf
	//	 	printf("random inf is %ld", random_inf) ;
		 	int counter4 = 0 ;
		 	while(current_animal!=NULL) // until it reaches till the end
		 	{		
		 	if(disease==1&&(current_animal->disease_status==2||current_animal->disease_status==1)) //if FMD infectious: for now if TB no passive detection occurs, if needed, add this: "(disease==0&&current_animal->disease_status>=2)||(" before hand and ")" in the end
			  {
			  
			  if(detection_scenario==0||detection_scenario==1)
			  {
			  	
			  			
			  	// Change the disease status and number recorded in FarmProductionStatus
			  						  	
					if(current_animal->disease_status==1)
					{
						FarmProductionStatus[pro_id][2] -- ;
						FarmProductionStatus[pro_id][4]++ ;
						FarmDataStatus[(int)floor(pro_id/3)][3]++;
						}	
			  		else if(current_animal->disease_status==2)
			  		{ //infectious
			  			FarmProductionStatus[pro_id][3]--; //FMD fectious minus 1
			  			FarmProductionStatus[pro_id][4]++ ; //One increase in the recovered
			  			FarmDataStatus[(int)floor(pro_id/3)][2]--;
			  			FarmDataStatus[(int)floor(pro_id/3)][3]++;
			  			long farm_id = (int)floor(pro_id/3) ;
			  			//printf("isolate is %d",FarmData[farm_id][5]);
			  		//	printf("id is %d",FarmData[(int)floor(pro_id/3)][5]);
			  		//system("pause");
			  			if(FarmDataStatus[farm_id][2]!=(FarmProductionStatus[farm_id*3][3]+FarmProductionStatus[farm_id*3+1][3]+FarmProductionStatus[farm_id*3+2][3]))
			  			{
			  			//	printf("L3569: FD is %lld FPS %f - %f - %f",FarmData[farm_id][7],FarmProductionStatus[farm_id*3][3],FarmProductionStatus[farm_id*3+1][3],FarmProductionStatus[farm_id*3+2][3]) ;
						  	system("pause") ;
						  }
			  			
			  			
			  			if(FarmDataStatus[(int)floor(pro_id/3)][2]<0)
			  			{
			  			//	printf("L3503:now it's negative %lld", FarmData[(int)floor(pro_id/3)][7]) ;
			  				system("pause") ;
						  }
			  			
			  			if(counter4==random_inf)
			  			{ //isolation of samples
			  			//	printf("err5")	;
			  			counter4++;
			  		    if(FarmDataStatus[farm_id][5]==0)
 			            {
 		                current_animal->ptr_isolate->isolated = 1 ; //isolated, ok what happens if this animal moved after this?
			  			current_animal->ptr_isolate->isolated_time = today_date ;
			  			(*nth_isolate_detected)++;
		               	FarmDataStatus[farm_id][5]=1;
						   if(rid_status[(int)FarmData[farm_id][5]]==0)
						   {
						   	rid_status[(int)FarmData[farm_id][5]]=1 ;
							   }	
			            }
			  			
			  			if(*nth_isolate_detected>=num_isolate_upper)
			  			{
			  				*sample_limit = 1 ;
			  				printf("stopped");
			  				return 0;
						  }
			  		
			  		//	printf("err6")	;
						}
						else
						{
							counter4++;
						}
					  }
			  			
			  		current_animal->disease_status= 3; // becomes immune - then have to cancel the event for them to beome immune or just skip it by saying only if no immune , immune occurs	
			  		//current_animal->ptr_isolate = NULL ;
			  	//	(*nth_isolate)++; // increment isolate counter by 1
				  }
			  }
			  current_animal = current_animal->next_node ; //next animal
			  } //while loop done
			  FarmDataStatus[farm_id][1] = 1 ; //becomes detected
			  	// and add a movement-ban-lift a month later
			  	if(floor(today_date)+movement_ban_days<sim_days)
			  	{
			  	struct event_node *adding_new_event;
			  			adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      			adding_new_event->akey=-1; // akey does not matter
		      			adding_new_event->src_pro_id = pro_id; 
		      			adding_new_event->des_pro_id = -100;
		      			adding_new_event->src_testarea=-100;
			  			adding_new_event->event_type=10;//movement ban lift
			            int day = floor(today_date)+movement_ban_days ;
			 			struct event_node *current_node1;
                        current_node1 = event_day[day];
                        if(current_node1 == NULL)
                         {
                         event_day[day] = adding_new_event;
						 }
                         else
                         {
                         adding_new_event -> next_node = current_node1;
                         event_day[day] = adding_new_event;
                         }	
				  }
			  	
		}
	//	printf("yep done") ; // else1 done
			  		
			/* Below is when we need to think about passive regular surveillance in addition to slaughterhouse inspection
				  	else if(disease==0)
				  	    {
				  	   	 if(current_animal->disease_status==2)
						    {
							FarmProductionStatus[pro_id][3]--; 
						    }
						 else if(current_animal->disease_status==3)
						    {
							FarmProductionStatus[pro_id][4]--; 
						    }
					     	FarmProductionStatus[pro_id][0]--; // total number -1 
						
						// remove this animal from the FarmProductionList
						previous_animal = current_animal->previous_node ; //get the previous animal
						next_animal = current_animal->next_node ; //get the next animal so that they can be rewired
						current_animal->previous_node = NULL ; //reset to NULL
						current_animal->next_node = NULL ; //reset to NULL
						
						previous_animal->next_node = next_animal ;
						next_animal->previous_node = previous_animal ;
						
						//then add this animal to slaughtered list
							if(inf_slaughtered[today_date]==NULL)
							{
								inf_slaughtered[today_date] = current_animal; //
							}
							else
							{
								next_slaughtered = inf_slaughtered[today_date] ;
								inf_slaughtered[today_date] = current_animal ;
								current_animal->next_node = next_slaughtered ;
							}
							current_animal = next_animal ; //next animal
						 }
			======= up to here should be restored when TB passive regular surveillance needs to be considered*/			 
    //  printf("detection done") ;  
    } // END OF DETECTION
	
	/*===================================DETECTION PART ENDS=========================================================================*/
	} //	if (next_non_markov_date>day_to_markov+today_date) DONE
	else
	{
		today_date = next_non_markov_date;
	}
//ASSESSING MARKOV OR NOT DONE	
return(today_date) ;
    
    
}//END OF FUNCTION
/*------------------------------------------------------------------------------------------------*/
/*COUNT NUMBER OF INFECTED AND DETECTED FARMS*/
/*--------------------------------------------------------------*/
void count_farms_infected_detected(double **FarmDataStatus, double **OutPut, long num_total_farms2, int sim)
{
long i;
OutPut[sim][0]=0;
OutPut[sim][1]=0;
for (i=0;i<num_total_farms;i++)
{
	if(FarmDataStatus[i][2]>0||FarmDataStatus[i][3]>0)
	{
	 OutPut[sim][0]++; //infected	
	 //printf("one more infected");
	}
	if(FarmDataStatus[i][1]==1)
	{
	 OutPut[sim][1]++;//detected
	// printf("one more detected");
	}
}
 	
	
	
}
/*----------------------------------------------------------------------------------------------------------*/
/* READ SEQUENCE DATA*/
/*------------------------------------------------------------------------------*/

void read_sequence_data(char SequenceData[], int seq_length, char seq_master[])
{
     FILE *Seq = fopen(SequenceData, "r");
     int c;
     int i = 0;
     if (Seq) {
      while ((c = fgetc(Seq)) != EOF)
                {
    
              // printf("c is %c\n", c);
              if (isspace(c)) // this was necessary because for some reasons the sequence text file contains empty character
                  {
              //     printf("this is empty");
                   }
                else{
                     seq_master[i]=c;
                 //  printf("seq is %c\n", seq_master[i]);
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


/*-----------------------------------------------------------------------------*/
/*Export CSV file of the number of bTB detected\infected data*/
/*------------------------------------------------------------------------------*/
int write_OutPut(char OutDataFile[],double **OutPut, int tot_iterations, int num_OutPut)
{

	FILE *Out = fopen(OutDataFile,"w");
	int line_num, col_num;
	
	for (line_num = 0 ; line_num < tot_iterations; line_num ++)
	{
		for (col_num = 0;col_num <num_OutPut ; col_num++ )
		
		{
	    fprintf(Out,"%lf,",OutPut[line_num][col_num]);
        }
        fprintf(Out,"\n");
    }
	fclose(Out);
	return 0;
}

/*----------------------------------------------------------------------------------------------*/
/*Export Sequence data---------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------*/
int export_IsolateData(char* OutIsolateDataFile, char* OutNexusFile,char* OutNexusFile_E,char* OutNexusFile_R,char* OutNexusFile_H_I,char* OutNexusFile_H,struct isolate** head_list_isolate, int **MembershipData, int num_max_isolate,int* nth_isolate_detected,int num_isolate_upper, int seq_length, int year, double** isolate_table,double **hid_rid_Data, int rid_inf_array[],int num_active_herds_rid_array[], int n_rid2)
{
	FILE *Out = fopen(OutIsolateDataFile, "w") ;
	FILE *Nexus = fopen(OutNexusFile, "w") ;
	int line_num, id,member ;
	int i = 0;
	int j,k,l,m ;
	int temp_count = 0;
	int isolate_count = 0;
	int rid = 0;
	double unique_rid_array[16] = {0} ;
	int rid_counter = 0 ;
	line_num = 0;
	double value1,value2 ;
	double min_t = 3000 ; 
	double max_t = 0 ;
	fprintf(Nexus,"#nexus\n");
	fprintf(Nexus,"begin data;\n");
	fprintf(Nexus,"Dimensions NTax=%d NChar =%d;\n",*nth_isolate_detected,seq_length-1);
	fprintf(Nexus,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Nexus,"matrix\n");
	
		while(head_list_isolate[line_num]!=NULL && line_num < num_max_isolate && isolate_count < num_isolate_upper)
		{
		fprintf(Out,"%d,",head_list_isolate[line_num]->id);
		fprintf(Out,"%d,",head_list_isolate[line_num]->parent_id);
		fprintf(Out,"%f,",head_list_isolate[line_num]->d_time);
		fprintf(Out,"%f,",head_list_isolate[line_num]->isolated_time);
		fprintf(Out,"%d,",head_list_isolate[line_num]->num_mut);
		fprintf(Out,"%lld,",head_list_isolate[line_num]->pro_id);
		fprintf(Out,"%d\n",head_list_isolate[line_num]->isolated);	
		
		if(head_list_isolate[line_num]->isolated==1)
		{
			i = 0;
			id = head_list_isolate[line_num]->id;
			member = MembershipData[(int)floor((head_list_isolate[line_num]->pro_id)/3)][year];
			rid = hid_rid_Data[(int)floor((head_list_isolate[line_num]->pro_id)/3)][1] ;
			if(unique_rid_array[rid-1]==0)
			{
				unique_rid_array[rid-1] = 1 ;
				rid_counter++;
			}
			
		value1 = roundf((head_list_isolate[line_num]->d_time)*100)/100 ;
		value2 = roundf((head_list_isolate[line_num]->isolated_time)*100)/100 ;
		
		if(value2>0)
		{
			fprintf(Nexus,"Isolate%d_%.2f_%.2f_%d ",id, value1,value2,member) ;
			isolate_table[isolate_count][0] = id ;
			isolate_table[isolate_count][1] = rid; // this is the original rid
			isolate_table[isolate_count][2] = value2 ;
			if(value2<min_t)
			{
				min_t = value2;
			}
			if(value2>max_t)
			{
				max_t = value2;
			}
			 
		while((head_list_isolate[line_num]->sequence[i])!='\0')
        {
        	fprintf(Nexus,"%c",head_list_isolate[line_num]->sequence[i]) ;
			i++ ;
		}
		fprintf(Nexus,"\n");
		}
		isolate_count++;
	
		
		}
		line_num ++ ;
		}
		
		
		
	//	int i = 0;
	//	for (col_num = 0;col_num <num_isolate_vars ; col_num++ )
		
	//	{
	 //   fprintf(Out,"%f,",IsolateData[line_num][col_num]);
     //   }
      //  while(IsolateSeqData[line_num][i]!='\0')
     //   {
      //  	fprintf(Out,"%c",IsolateSeqData[line_num][i]) ;
	//		i++ ;
		
        
      //  fprintf(Out,"\n");
    
    fprintf(Nexus,";\n");
	fprintf(Nexus,"end;");
	fclose(Nexus);
    fclose(Out) ;
    
   // printf("Y4") ;
/*From here open Nexus file for each sub-sampling and export them*/
    int interval_n = floor(isolate_count/rid_counter)+1 ;
    double interval_t = (max_t - min_t)/interval_n ;
    int sum_inf_array[16] = {0} ;
    int sum_active_array[16] = {0} ;
    int id_array[150] = {-1} ;
    int random_int ;
    int count_E = 0;
    FILE *Nexus_E = fopen(OutNexusFile_E, "w") ;
    FILE *Nexus_R = fopen(OutNexusFile_R, "w") ;
    FILE *Nexus_H_I = fopen(OutNexusFile_H_I, "w") ;
    FILE *Nexus_H = fopen(OutNexusFile_H, "w") ;
    
   
    for(i=1 ; i <= 16; i++)
    {
    	if(unique_rid_array[i-1]==1) //means there is a sample
    	{
    	for(k = 0; k < interval_n ; k++)
    	{
    		temp_count = 0;
    		// choose id that is in a given time stamp and in this region, how to pick up and save id?
    		for(j=0; j < isolate_count; j++)
    		{
    			if((isolate_table[j][1]==i&& k< interval_n - 1 && (isolate_table[j][2]>=min_t+k*interval_t && isolate_table[j][2]<min_t+(k+1)*interval_t))||(isolate_table[j][1]==i&&k==interval_n - 1 && (isolate_table[j][2]>=min_t+k*interval_t && isolate_table[j][2]<=min_t+(k+1)*interval_t)))
				{
					isolate_table[j][3] = 1 ;
					id_array[temp_count] = j ;
					temp_count++ ;
				 } 
			
			}
			//then choose one sample
			if(temp_count>0)
			{
				random_int = rand()%temp_count ;
				isolate_table[id_array[random_int]][3] = 2 ;
				count_E++;
				//immediately store the data? No do it later.
				// so only need to choose those have 2 in column 3
			}
		}
	
		}
	} //sampling E done
//	printf("Y5") ;
	//exporting Nexus_E
	fprintf(Nexus_E,"#nexus\n");
	fprintf(Nexus_E,"begin data;\n");
	fprintf(Nexus_E,"Dimensions NTax=%d NChar =%d;\n",count_E,seq_length-1);
	fprintf(Nexus_E,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Nexus_E,"matrix\n");
	
	//sampling R starts
	i = 0;
	while(i<count_E)
	{
		random_int = rand()%isolate_count ;
		if(isolate_table[random_int][4]==0)
		{
		isolate_table[random_int][4]= 2 ;
		i++ ;	
		}
		
		
	} //sampling R done
//	printf("Y6") ;
	fprintf(Nexus_R,"#nexus\n");
	fprintf(Nexus_R,"begin data;\n");
	fprintf(Nexus_R,"Dimensions NTax=%d NChar =%d;\n",count_E,seq_length-1);
	fprintf(Nexus_R,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Nexus_R,"matrix\n");
	
	int sum_I = 0 ;
	int sum_H = 0;
	int random_I ;
	int random_H ;
	int rand_temp ;
	//sampling H_I starts
	// need the number of detected herds in each region
	for(i=0;i<n_rid2;i++)
	{
		sum_I = sum_I + rid_inf_array[i] ;
		sum_inf_array[i] = sum_I ;
		
		sum_H = sum_H + num_active_herds_rid_array[i];
		sum_active_array[i] = sum_H ;
		
	}
//	printf("Y6-2") ;
	// now choose samples for the number that is required
	i = 0;
	while(i<count_E)
	{
	  random_I = rand()%sum_I ;
    j = 0;
    while(random_I>=sum_inf_array[j])
    {
    	j++;
	}
	//then choose from region J
		if(rid_inf_array[j]>0)
	{
		rand_temp = rand()%rid_inf_array[j];
	//	printf("rand is %d", rand_temp) ;
		k = 0;
		l = -1;
		while(k<rand_temp)
		{
			l++;
	     if(isolate_table[l][1]==j+1) //j is rid -1
			{
				k++;
			}
          
		}
		if(l==-1)
		{
			l=0 ;
		}
		if(isolate_table[l][5] == 0)
		{
			isolate_table[l][5] = 2 ;
			i++ ; //counter increment and sample taken
		}
		
	} //sampling H_I done
     }
//	printf("Y7") ;
	fprintf(Nexus_H_I,"#nexus\n");
	fprintf(Nexus_H_I,"begin data;\n");
	fprintf(Nexus_H_I,"Dimensions NTax=%d NChar =%d;\n",count_E,seq_length-1);
	fprintf(Nexus_H_I,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Nexus_H_I,"matrix\n");
  
  //sampling H starts
  
  	i = 0;
	while(i<count_E)
	{
	//	printf("i is %d", i) ;
	  random_H = rand()%sum_H ;
    j = 0;
    while(random_H>=sum_active_array[j]) //this is choosing region proportional to the number of active herds
    {
    	j++;
	}
	//then choose from region J
	if(rid_inf_array[j]>0)
	{
		rand_temp = rand()%rid_inf_array[j]; //right if it's choosing region that does not have any isolates...
	//	printf("rand is %d", rand_temp) ;
		k = 0;
		l = -1;
		while(k<rand_temp)
		{
			l++;
	     if(isolate_table[l][1]==j+1)
			{
				k++;
			}
          
		}
		if(l==-1)
		{
			l=0;
		}
		if(isolate_table[l][6] == 0)
		{
			isolate_table[l][6] = 2 ;
			i++ ; //counter increment and sample taken
		}
	}

		
	} //sampling H done
//	printf("Y8") ;
  //export data
 	fprintf(Nexus_H,"#nexus\n");
	fprintf(Nexus_H,"begin data;\n");
	fprintf(Nexus_H,"Dimensions NTax=%d NChar =%d;\n",count_E,seq_length-1);
	fprintf(Nexus_H,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Nexus_H,"matrix\n");
    int sample_id = 0;
    //go through isolate_table and write data in
    for(i=0;i<isolate_count;i++)
    {
    	
    	if(isolate_table[i][3]==2) //E
    	{
    		j = 0;
    		sample_id = (int)isolate_table[i][0] ;
    	fprintf(Nexus_E,"Isolate%d_%.2f_ ",sample_id, isolate_table[i][2]) ;	
    	while((head_list_isolate[sample_id]->sequence[j])!='\0')
        {
        	fprintf(Nexus_E,"%c",head_list_isolate[sample_id]->sequence[j]) ;
			j++ ;
		}
		fprintf(Nexus_E,"\n");
		}
		if(isolate_table[i][4]==2) //R
    	{
    			j = 0;
    		sample_id = (int)isolate_table[i][0] ;
    	fprintf(Nexus_R,"Isolate%d_%.2f_ ",sample_id, isolate_table[i][2]) ;	
    	while((head_list_isolate[sample_id]->sequence[j])!='\0')
        {
        	fprintf(Nexus_R,"%c",head_list_isolate[sample_id]->sequence[j]) ;
			j++ ;
		}
		fprintf(Nexus_R,"\n");
		}
		
		if(isolate_table[i][5]==2) //H_I
    	{
    			j = 0;
    		sample_id = (int)isolate_table[i][0] ;
    	fprintf(Nexus_H_I,"Isolate%d_%.2f_ ",sample_id, isolate_table[i][2]) ;	
    	while((head_list_isolate[sample_id]->sequence[j])!='\0')
        {
        	fprintf(Nexus_H_I,"%c",head_list_isolate[sample_id]->sequence[j]) ;
			j++ ;
		}
		fprintf(Nexus_H_I,"\n");
		}
		
		if(isolate_table[i][6]==2) //R
    	{
    			j = 0;
    		sample_id = (int)isolate_table[i][0] ;
    	fprintf(Nexus_H,"Isolate%d_%.2f_ ",sample_id, isolate_table[i][2]) ;	
    	while((head_list_isolate[sample_id]->sequence[j])!='\0')
        {
        	fprintf(Nexus_H,"%c",head_list_isolate[sample_id]->sequence[j]) ;
			j++ ;
		}
		fprintf(Nexus_H,"\n");
		}
    	
	}
	fprintf(Nexus_E,";\n");
	fprintf(Nexus_E,"end;");
	fclose(Nexus_E);
	fprintf(Nexus_R,";\n");
	fprintf(Nexus_R,"end;");
	fclose(Nexus_R)	;
	fprintf(Nexus_H_I,";\n");
	fprintf(Nexus_H_I,"end;");
	fclose(Nexus_H_I)	;
	fprintf(Nexus_H,";\n");
	fprintf(Nexus_H,"end;");
	fclose(Nexus_H)	;
   // printf("Y9") ;
    return 0 ;
	
}

/*-------------------------------------------------------------------------------------
/* EXPORT SEQUENCE IN NEXUS FORMAT------------------------------------------------------
/*-------------------------------------------------------------------------------------*/
/*
int export_seq_nexus(char* OutNexusFile, int* nth_isolate_detected, int seq_length,double **IsolateData, char *IsolateSeqData[])
{
	int line_num, col_num,i,id;
	double value1,value2 ;
	FILE *Out = fopen(OutNexusFile,"w") ;
	fprintf(Out,"#nexus\n");
	fprintf(Out,"begin data;\n");
	fprintf(Out,"Dimensions NTax=%d NChar =%d;\n",*nth_isolate_detected,seq_length-1);
	fprintf(Out,"Format DataType=DNA missing=? gap=-;\n");
	fprintf(Out,"matrix\n");
	for(line_num=0; line_num < *nth_isolate; line_num++)
	{
		i = 0 ;
		if(IsolateData[line_num][7]==1)
		{
		//if(IsolateData[line_num][7]==1 //only if isolated
		//{}
		id = roundf(IsolateData[line_num][4]) ;
		value1 = roundf(IsolateData[line_num][3]*100)/100 ;
		value2 = roundf(IsolateData[line_num][2]*100)/100 ;
		if(value2>0)
		{
			fprintf(Out,"Isolate%d_%.2f_%.2f ",id, value1,value2) ;
		while(IsolateSeqData[line_num][i]!='\0')
        {
        	fprintf(Out,"%c",IsolateSeqData[line_num][i]) ;
			i++ ;
		}
		fprintf(Out,"\n");
		}	
		}
		
			
	}
	fprintf(Out,";\n");
	fprintf(Out,"end;");
	fclose(Out);
	return 0;
	
}
*/
