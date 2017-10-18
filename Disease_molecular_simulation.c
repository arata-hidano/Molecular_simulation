/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEST VERSION 1. 
1. USE SMALL DATA TO CHECK IF MUTATION AND TRANSMISSION ARE OCCURRING.
2. NO BIRTH, NO SURVEILLANCE TESTING.
3. A SHORTER SIMULATION.
Check //@@@ part because these parameters were changed
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

/*DEFINE VARIABLES*/
#define num_production_type 3
#define num_total_farms 45965  // this is the total number of herds that appeared in movement\longvity table
int num_total_farms2 = num_total_farms ;
int const num_animals = 3624420;  // as of 2000 July 1st
int const num_total_animals = 16534951;  
#define sim_years 1 //@@@ 
int num_max_isolate = 200 ;  // should change this
int seq_length_FMD = 633 ;
#define seq_length_TB 700 //@@@length of TB gene to be checked

/*VARIABLES TO CHANGE PARAMETERS*/
int seq_length = 634;
/*MUTATION MODEL OR NOT?*/
int mutation_model = 1; //1 if mutation model
int detection_scenario = 0 ;
/*WHICH DISEASE?*/
int disease = 1; //0 if FMD, 1 if FMD
long long num_farm_production = num_total_farms*num_production_type ; // This is the ID for each farm production unit
int num_farm_var = 11 ; // 11 variables for each farm; [0] farm_id, [1] x, [2] y, [3] testarea_num, [4] disease_status, [5] testing schedule, [6] tb_detected, [7] sum of occult(reactive) for TB, infectious for FMD, [8] sum of infectious for TB, immune for FMD
                        // [9] farm id that infected this farm, [10] date that this source farm infected this destination farm

char FarmDataFile[] = "all_herds_details_export.csv";
int column_prostatus = 9 ; // number of variable in the table of infection status
long long current_pro_id ;
// at this monemnt; [0] N, [1] sus, [2] exposed, [3] occult, [4] detectable,[5] infection pressure, [6] cumu_pressure, [7] mutation_pressure, [8] move_ban

/* Variables related to individual animal*/

int num_animals_var = 8; // [0] akey, [1] farm_id, [2] age_day, [3] age_type, [4] sex, [5] breed, [6] pregnant_status, [7] disease_status
// just need to initialise [7] and [8]
char AnimalDataFile[] = "tanimals_present_1july2000_toy.csv";

/* Variables related to events(movement, birth,test)*/
char MoveDataFile[] = "tmovements_since2000july_export2.csv"; ///need to make it short
long num_moves = 9681150;
//long num_moves = 888677;
int num_moves_vars = 6; // serial_akey, date, src_serial_hid, des_serial_hid, age_type, src_testarea

/* //@@@ IGNORE BIRTH
char BirthDataFile[] = "/C_run/birth_table_since2000july.csv";*/
int num_births = 9667100 ;
int num_births_vars = 5;// akey, bdate, src_farm, sex, breed 

/* //@@@ IGNORE TEST
char TestDataFile[] = "/C_run/tb_test_schedule_export.csv";*/
int num_tests = 742 ;
int num_tests_vars = 2;


/* Variables related to simulations*/

int sim_days= 365*sim_years; // X years
//341 farms detected TB in this 3 years
long long i, j;
int iteration;
int tot_iterations = 1;
int num_OutPut = 2; //infected and detected
char OutPutFile[] = "OutPut2.csv";
/* Variables related to disease*/
int day_S3toS4, day_S2toS3;
int max_S3toS4_TB = 0.4*365; // max day to detect - 1 Conlan(2014) modified
int max_S2toS3_TB = 5*365; // max day to occult - 1  Brooks-Pollock's parameter 11.1
int max_S3toS4_FMD = 20 ;
int max_S2toS3_FMD = 5 ;
int max_S2toS3, max_S3toS4, day_to_S3, day_to_S4 ;
double Se_occult = 0.2;
double Se_detect = 0.45;
double Se_slaughter = 0.3 ;// @@@@ have to revise this value
double beta_a = 1;//within-herd transmission parameter
double beta_b = 0.00015;//wildlife transmission parameter
double mu = 2.7*0.001*633/365; // use substitution rate used in Hall, which was further converted to the substitution rate for the whole VP1 seq and per day
 //must be random draw from some distribution 0.003 - 0.036 Brooks-Pollock
int *current_id_isolate; // a pointer to a memory that keeps the track of isolate id
char Seq_master[634]; //Master (initial) sequence for FMD VP1

/*File names*/
char SequenceData[] = "FJ785304_1_VP1.txt"; //Seq file for FMD
char OutputXML[] = "test.xml"; //Output XML file name
const char* OutIsolateDataFile_origin = "Sequence_data_iteration";
const char* extension = ".csv" ; 

int num_isolate_vars = 8;
int num_isolate_seq_vars = 2;
// [0]: akey, [1]: pro_id, [2]: isolated_time, [3]: diverged_time, [4]: id, [5]: parent_id, [6]: num_mut, [7]: isolated (1==Yes).
int *nth_isolate;
int *nth_slaughter;

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


//seq_length = seq_length_FMD;
struct isolate {
		char sequence[634];
		int id;
		int parent_id ; //id of parent
		double d_time ; //day this isolate diverged from the parent
		int num_mut ; // number of mutations
					};
   
/* Array that points to each detected isolates*/

struct isolate *head_list_isolate[200]; 

/* FUNCTION DEFINITIONS */
  
void read_farm_data(); //Function for reading the data
void read_movement_data() ;
void read_animal_data();
void read_test_data();
void add_event_node();
void add_animal_node();
void visualize_list();
void visualize_animals() ;
void move_animal_unit();
void read_birth_data();
void test_farms();
double update_markov_date();
void count_farms_infected_detected();
int write_OutPut();
void move_production_type();
 
void read_sequence_data() ;
int createdata();
void JukesCantor();

/* Variables related to farm and FarmProductionStatus*/
//int num_farms = 16950; // number of farms present at 2000 july 




/*-------------------------------------------------------------------------------*/


/* what's the problem about not tracing non-infected animals?
For instance, we can update the status of a farm when they receive bTB infected animals
and then adjust the age structure and herd size of farm.
By doing this we assume herd size stays the still and age structure is the same.
But can we really know if the herd size is going to be increased or not?
If not this assumption seems to be ok*/

/* MAIN PROGRAM */
int main(void){
	
       
    
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
   	  	 long long **FarmData = (long long**)malloc( sizeof(long long *) * num_total_farms);
   	  	 for(i = 0; i < num_total_farms; i++)
      	  {
     	       /*Then, for each row, specifies the required memory size for column*/
     	     FarmData[i] = (long long*)malloc( sizeof(long long) * num_farm_var); 
			  /* what kind of farm variable do we need?*/ 
     	   }
      	
     	read_farm_data(FarmDataFile, FarmData, num_total_farms);
		 //printf("farm data read");     
/*=====2.3 Read in animal Data===========================================================*/
       		long long **AnimalData = (long long**)malloc( sizeof(long long *) * num_animals);
   	  		 for(i = 0; i < num_animals; i++)
      	 	 {
     	       /*Then, for each row, specifies the required memory size for column*/
     	     AnimalData[i] = (long long*)malloc( sizeof(long long) * num_animals_var);      	   		
			}
					
      	read_animal_data(AnimalDataFile,AnimalData,num_animals) ;
      	printf("The first animals is %lld",AnimalData[0][0]) ;
      	//printf("animal data read");
      /* Check if it reads AnimalData properly
	  	printf("Animal id is %lld",(long long) AnimalData[0][0]); // intersting it's now working, why do we need declaration?
        long long temp = AnimalData[0][0];
        int temp2 = AnimalData[0][1];
        printf("Animal id and farm is %lld, %d", temp, temp2);
      	system("pause") ;*/

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
   	double **TestData = (double**)malloc( sizeof(double *) * num_tests);
   	    	for(i = 0; i < num_tests; i++)
   	    	 {
   	    	 TestData[i] = (double*)malloc( sizeof(double) * num_tests_vars);  
   		     	}
   	    //@@@	read_test_data(TestDataFile, TestData, num_tests);
   	    //	printf("test data read");
   	    	
   	    	
/*2.6 Read birth data*/
    long long **BirthData = (long long**)malloc( sizeof(long long *) * num_births);
   	    	for(i = 0; i < num_births; i++)
   	    	 {
   	    	 BirthData[i] = (long long*)malloc( sizeof(long long) * num_births_vars);  
   		     	}
   	   //@@@      read_birth_data(BirthDataFile, BirthData, num_births) ;
   	      // printf("birth data read");
		   
/*2.7 Read sequence data*/
read_sequence_data(SequenceData,seq_length_FMD,Seq_master) ;
 	         

/*Create IsolateData frame that stores isolate information*/

double **IsolateData = (double**)malloc( sizeof(double *) * num_max_isolate);
char* IsolateSeqData[num_max_isolate] ; //= (char*)malloc(sizeof(char*)*num_max_isolate);
   	    	for(i = 0; i < num_max_isolate; i++)
   	    	 {
   	    	 IsolateData[i] = (double*)malloc( sizeof(double) * num_isolate_vars);  
   	    	 IsolateSeqData[i] = (char*)malloc(sizeof(char)*634) ;
   		     	}
	//printf("IsolateData created");

/*Create OutPut data*/
    double **OutPut = (double**)malloc( sizeof(double *) *tot_iterations);
   	    	for(i = 0; i < tot_iterations; i++)
   	    	 {
   	    	 OutPut[i] = (double*)malloc( sizeof(double) *num_OutPut); 
				   		     }
				   		     
/*Create FarmProductionStatus*/
double **FarmProductionStatus = (double**)malloc( sizeof(double *) * num_farm_production);
for(i = 0; i < num_farm_production; i++)
      	  	{
     	       /*Then, for each row, specifies the required memory size for column*/
     	     FarmProductionStatus[i] = (double*)malloc( sizeof(double) * column_prostatus); 
     	    }
      	/*============================================================================================*/
 
/* Set memory for linked lists*/      
struct animal_node **animal_node_pointer = (struct animal_node**)malloc( sizeof(struct animal_node*) * num_total_animals);	;
for(i=0; i < num_total_animals; i++)
{
	animal_node_pointer[i] = (struct animal_node*)malloc(sizeof(struct animal_node)) ;
	//animal_node_pointer[i]= NULL ;
	animal_node_pointer[i]->akey = -1; // initialisation
}
//animal_node_pointer = (struct animal_node*)malloc( sizeof(struct animal_node) * num_total_animals);	
// animal_node_pointer points to the address of the first animal of the array with a length of num_total_animals
//struct animal_node *animal_node_pointer = (struct animal_node*)malloc( sizeof(struct animal_node) * num_total_animals);
struct animal_node *inf_slaughtered[sim_days] ; // a pointer to animal_node that stores slaughtered and infected animals
// maybe need to add slaughtered and non-infected linked list too

struct animal_node* FarmProductionList[num_farm_production]; // pointer to first animal at each farm-production
struct inf_move *Farm_ptr_inf_move[num_total_farms] = {NULL} ; // this records moves that had infected animals from the source farm
struct event_node* event_day[sim_days]; // a pointer to a struct of event

 // if(Farm_ptr_inf_move[37761]==NULL)
 // {
  //	printf("Say yes it's null'") ;
  //}
  //else
 // {
  //	printf("37761 is not null") ;
  //}
  
  /*Parameters*/
  if(disease==0) //TB
      	{
      		max_S2toS3 = max_S2toS3_TB;
      		max_S3toS4 = max_S3toS4_TB ;
		}
		else if(disease==1)
		{
			max_S2toS3 = max_S2toS3_FMD ;
			max_S3toS4 = max_S3toS4_FMD ;
		  } 
		  
	 
    
    
    
/*===============================================================================================================================================================================================================
------------------START OF THE SIMULATION-----------------------------------------------------------
============================================================================================*/  



/*===========START OF ITERATION=======================================================================*/
for(iteration=0; iteration<tot_iterations; iteration++)
{//iteration starts
		    	
      printf("Iteration started");
      long int random_id = rand()%num_animals;
      int var1 = 0 ;
      nth_isolate = &var1;
      
      int var2 = 0;
      nth_slaughter = &var2;
      printf("nth_slaughter is %d", *nth_slaughter);
      char* OutIsolateDataFile;
      
      char number[1] ;
      
      snprintf(number, 1, "%d", iteration) ;
      OutIsolateDataFile = malloc(strlen(OutIsolateDataFile_origin)+1+1+4) ; //size of _origin, iteration, /0(null terminate), and extension(csv)
      strcpy(OutIsolateDataFile, OutIsolateDataFile_origin);
      
      strcat(OutIsolateDataFile, number) ;
      strcat(OutIsolateDataFile, extension) ;
	  srand((unsigned)time(NULL));	
       
      	
      	int var3 = 0;
		current_id_isolate = &var3; //reset the counter for current id for isolate
      	// INITIALISE THE FARM TB DETECTED STATUS, SUM OF OCCULT, SUM OF DETECTABLE
      	for(i=0; i< num_total_farms; i++)
      	{
      		FarmData[i][6] = 0;
      		FarmData[i][7] = 0;
      		FarmData[i][8] = 0;
			FarmData[i][9] = -1;
			FarmData[i][10] = -1; // THESE are necessary because they need to be initialised after each iteration
      		// initialise Farm_ptr_inf_move 
      		//Farm_ptr_inf_move[i]->des_farm_id = -1 ;
		  }
		  
		  for(i = 0; i < num_max_isolate; i++)
   	    	 {
   	    	 	for (j=0; j < num_isolate_vars; j++)
				{
				IsolateData[i][j] = 0 ; 
				}
				
				IsolateSeqData[i] = "";
				
   	    	  
   		    }
   		     	
   		      for(i=0;i<tot_iterations;i++)
 			 {
  				for (j=0; j < num_OutPut; j++)
				{
				OutPut[i][j] = 0 ; 
				}
 			 }
      	
/* Following codes to make linked lists need to be repeated in each iteration*/
/*=====Prepare Vector of day "event_day" that links to EventNode that stores information for any events (move/born/culling, infection status change)*/

/*=====2.2.1 Farm-production pointer 1: Pointer to table of infection status=======================*/
	  	// This is table storing number of infected, susceptible and so on for each management unit
	  	
	  		printf("4");
   	  	 	for(i = 0; i < num_farm_production; i++)
      	  	{
     	       /*Then, for each row, specifies the required memory size for column*/
     	    // FarmProductionStatus[i] = (double*)malloc( sizeof(double) * column_prostatus); 
     	     FarmProductionList[i] = NULL; // initialise the animal struct
			  /* what kind of farm variable do we need?*/ 
			  for (j=0; j < column_prostatus; j++)
			  {
			  	FarmProductionStatus[i][j] = 0 ; 
				  // initialise the status table, but keep herd id 
			  }
			  // then read the initial bTB status and total N
     	   	}
     	
     	/*====2.2.2 Farm-production pointer 2: Pointer to animal array===========================*/   	 		  
          for(i = 0; i < sim_days; i++)
                {
                event_day[i] = NULL;
                inf_slaughtered[i] = NULL ;
                }
                printf("5") ;
      //	printf("event is fine");
      	/* Create a vector of size of number of animals that store a pointer to animal_node*/
      	
      	
		  //for (i=0; i < num_total_animals; i++)
      //	{
      	//	animal_node_pointer[i] = (struct animal_node*)malloc( sizeof(struct animal_node)); 
      		 //initialise
      		//Farm_ptr_inf_move[i] = NULL ;
		 // }
       printf("animal_node_pointer is fine"); 
        
/*======================ADDING ANIMAL NODE=============================================================*/ 
          // Also here we can choose which animal is infected
          // Also we can count number of animas in each status if necessary to pass them to FarmProductionStatus
		  
		  //struct animal_node *new_node;
		  struct event_node *new_event;
		  // prepare a pointer that points to animal_node and event_node 
for (i=0; i < num_animals; i++)
          { 
                
                //new_node = (struct animal_node*)malloc(sizeof( struct animal_node )); 
                long long current_akey = (long long)AnimalData[i][0];
                //printf("current akey is %lld",current_akey); 
                animal_node_pointer[current_akey]->akey = current_akey;
                //printf("current akey is %lld",animal_node_pointer[current_akey]->akey);
                //new_node -> akey = current_akey ; // extract from animal data
                animal_node_pointer[current_akey]->age_day = (int)AnimalData[i][2];
                //new_node -> age_day = (int)AnimalData[i][2] ;
              //  if(current_akey==2454078)
               // {
               // 	printf("hi i'm here at %lld",(long long)AnimalData[i][1]);
			//	}
                
                long long current_farm = (long long)AnimalData[i][1] ; // farm id
               // int current_type = (int)AnimalData[i][3] ; // which age group (production unit)
                current_pro_id = AnimalData[i][1]*3 + AnimalData[i][3] ; // initialise the current_pro_id
                //current_pro_id = (int)AnimalData[i][8];
                animal_node_pointer[current_akey]->current_pro_id = current_pro_id;
				//new_node-> current_pro_id = (int)AnimalData[i][8]; //initial farm
                
                animal_node_pointer[current_akey]->sex = (int)AnimalData[i][4] ;
                //new_node -> sex = (int)AnimalData[i][4] ;
                animal_node_pointer[current_akey]->breed = (int)AnimalData[i][5];
                //new_node -> breed = (int)AnimalData[i][5];
            	animal_node_pointer[current_akey]->pregnant_status = (int)AnimalData[i][6];
				//new_node -> pregnant_status = (int)AnimalData[i][6];        
				
            	
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
        
        if (i == random_id)
        {
        	printf("This animal is infected") ;
        	system("pause") ;
        	animal_node_pointer[current_akey]->disease_status = 1;
        	//new_node -> disease_status = 1 ; // infected
        	struct isolate *new_isolate ;
        	new_isolate = (struct isolate*)malloc(sizeof( struct isolate )); //new_isolate is a pointer, so it's compatible to pass this to ptr_isolate
            new_isolate -> id = 0 ; // this is the seed
            strcpy(new_isolate->sequence,Seq_master);
            new_isolate -> parent_id = -1 ;
            new_isolate -> d_time = 0 ;
            new_isolate -> num_mut = 0 ;
            
            animal_node_pointer[current_akey]->ptr_isolate = new_isolate; 
        	//new_node -> ptr_isolate = new_isolate; //assigning a pointer new_isolate to a pointer ptr_isolate
        	FarmProductionStatus[current_pro_id][2]++; // increment the counter for infected
        	printf("%f animal is infected",FarmProductionStatus[current_pro_id][2] ) ;
        	system("pause") ;
			FarmData[current_farm][4] = 1; // farm is infected
        	FarmData[current_farm][9] = -1 ; // this is seed farm
        	FarmData[current_farm][10] = 0 ; // infected time is 0
        	*current_id_isolate = *current_id_isolate + 1;
		}
		else
		{
			animal_node_pointer[current_akey]->disease_status = 0 ;
			//new_node -> disease_status = 0 ;
			animal_node_pointer[current_akey]->ptr_isolate = NULL ;
			//new_node -> ptr_isolate = NULL ;
			FarmProductionStatus[current_pro_id][1]++;
		}
                
                  /*Can ADD OTHER INFECTED ANIMALS*/
          //set probability
          /*
          double p_inf = 0.0001;
          double random_value = ((double)rand()/(double)RAND_MAX) ;
          if (random_value < p_inf)
          {
          	double random_value = ((double)rand()/(double)RAND_MAX) ;
          	if(random_value<0.9)
          	{
          	new_node -> disease_status = 3;
            }
            else
            {
            new_node -> disease_status = 2;
			}
		  }
		  // 0: susceptible 1: exposed 2: occult 3: infectious
		  */
          /*This part is irrelevant unless occult condition is considered*/    
	// Exposed starts  
            // Now if the cattle is exposed status, choose time to be occult
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
            	printf("add S3");
            	system("pause") ;
                }
				/* Following was to check the bug in adding new event, which occurred because
            	I forgot to specify the size of new_event and then old memory was used
            	printf("My akey is %lld", current_akey) ;
            	system("pause") ;
            	printf("Day is %d", day_to_occult) ;
            	system("pause") ;
				visualize_list(event_day, day_to_occult);
				system("pause") ;*/
			}    //exposed up to here
            
			
			// If animal is in occult then choose a date to be infectious
			/* Occult starts
			if (new_node -> disease_status==2)
			{
				//randomly select from some distributions
				// for now draw a random value from uniform distribution
			
				day_to_detect = rand()%max_day_detect+1;
				FarmProductionStatus[current_pro_id][3]++; // increment by 1
				FarmData[current_farm][7] ++; //increase number of occult
				
				if(day_to_detect<sim_days)
				{
				
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				
				new_event -> akey = current_akey;
				new_event -> src_pro_id = current_pro_id ;
				new_event -> des_pro_id = -100 ;
				new_event -> event_type = 7; // 4 is occult to detectable
            	new_event -> next_node = NULL;
            	new_event -> src_testarea = -100 ;
            //	printf("add detect");
            	add_event_node(event_day,day_to_detect, new_event) ;
			    }
			}
            occult up to here*/
			// If animal is detectable
			/* Detectable starts here
			if (new_node -> disease_status==3)
			{
				
				FarmProductionStatus[current_pro_id][4]++; // increment by 1
				FarmData[current_farm][8]++;
				//printf("this farm has detectable");
				} 
				detectable up to here*/
			// if animal is susceptible
			/*susceptible starts here
			else if (new_node -> disease_status==0)
			{
		

		  
				FarmProductionStatus[current_pro_id][1]++; // increment by 1
		  
				} 
				susceptible ends here*/
			// finally increment for the total number of animals in this type
			FarmProductionStatus[current_pro_id][0] = FarmProductionStatus[current_pro_id][0]+1;
			
	/*-------------DISEASE STATUS BLOCK END--------------------------------------------------------------------*/
			
			// then initialise the next_node
			animal_node_pointer[current_akey]->next_node = NULL;
             //   new_node -> next_node = NULL; 
			 animal_node_pointer[current_akey]->previous_node = NULL ;  
                //new_node -> previous_node = NULL ;
            /* ADD THE NEW NODE TO THE ARRAY */
                 add_animal_node(FarmProductionList, current_pro_id, animal_node_pointer[current_akey] ) ;
                 
            /* ADD address of this animal to the pointer list*/
            //animal_node_pointer[current_akey] = new_node ;
           } 
           
           
          // printf("next is %lld previous is %lld",animal_node_pointer[2454078]->previous_node->akey,animal_node_pointer[2454078]->next_node->akey);
         // if(animal_node_pointer[2454078]->previous_node==NULL)
         // {
         // 	printf("YES NULL");
		 // }
		 // if(animal_node_pointer[0]->previous_node==NULL)
         // {
         // 	printf("YES NULL2");
		  //}
		  printf("adding animal done");
           //system("pause");
           
           /* Check if animals were added properly
           for (i=0; i < 2; i++)
           {
           	visualize_animals(FarmProductionList,i) ;
		   }
		   system("pause");
		   */
           
        /*    for (i = 0 ; i < 40; i++)
      {
          visualize_list(event_day, i);
          }
          system("pause");*/
      
/*======ADDING ANIMALS DONE=====================================================================================================*/      
      
     

   	    	// Movements are going to be probabilistic rather than explicitly 
   	    	// moving recorded animals?
   	    	// If so, we only need number of animals moved in each production type
   	    	
   	    //	printf("Movement data is like this %lld, %lld, %lld", (long long)MoveData[0][0], (long long)MoveData[1][0],(long long)MoveData[2][0]) ;
   	
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
	/* Add movements to event array*/
	 // each farm struct is also a pointer to a struct   
          //printf("MoveData akey is %lld, %lld, %lld, %lld", (long long)MoveData[0][0],(long long)MoveData[1][0],(long long)MoveData[2][0],(long long)MoveData[3][0]);
		  //system("pause");

/* Steps:
1. ADD MOVEMENT TO EVENT NODE/
2. ADD BIRTH TO EVENT NODE.

*/


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
			//	if(current_akey==2454078)
			//	{printf("hi day %d",current_day);
			//	}
                //printf("%d current akey is %lld",i, (long long)MoveData[i][0]);
                //system("pause");
                //printf("%d current akey is %lld",i, current_akey);
                //system("pause");
                long long current_pro_id = (long long)MoveData[i][2]*3 + current_event_type;
                
                new_event -> akey = current_akey;
                // printf("current akey is %lld", new_event -> akey);
                 //system("pause");
                new_event -> src_pro_id = current_pro_id ;
				new_event -> des_pro_id =  (long long)MoveData[i][3]*3+current_event_type; // destination farm
                new_event -> src_testarea = (int)MoveData[i][5] ;
				new_event -> next_node = NULL;   
         
                
                /* ADD THE NEW NODE TO THE ARRAY */
                 add_event_node(event_day, current_day, new_event) ;
                 if (animal_node_pointer[current_akey]->akey == -1)//if this animal does not exist initially
                 {
                 //struct animal_node *new_animal;
                 //new_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
                 animal_node_pointer[current_akey]->akey=current_akey;
              //   if(current_akey == 2751445)
              //   {
              //   	printf("ok akey is %lld", animal_node_pointer[current_akey]->akey) ;
              //   	system("pause") ;
			//	 }
                 //new_animal -> current_pro_id = current_pro_id ; //only first movement
                 animal_node_pointer[current_akey]->current_pro_id = current_pro_id ;
                 //new_animal->breed=6;//now just assume it's unknown
                 animal_node_pointer[current_akey]->breed = 6 ;
                 //new_animal->sex=1; //now assume all female but if needs to be precise, have to get the data from tlogevity and add to movement data
                 animal_node_pointer[current_akey]->sex=1;
				// new_animal->disease_status=0; //susceptible
				animal_node_pointer[current_akey]->disease_status=0;
				 //new_animal->num_births=0; //unknown
				 animal_node_pointer[current_akey]->num_births=0;
				 //ew_animal->type=current_event_type;
				 animal_node_pointer[current_akey]->type=current_event_type;
				 //new_animal->pregnant_status=0;
				 animal_node_pointer[current_akey]->pregnant_status=0;
				 //new_animal->ptr_isolate = NULL;
				 animal_node_pointer[current_akey]->ptr_isolate = NULL;
				 //new_animal->next_node=NULL;
				 animal_node_pointer[current_akey]->next_node=NULL;
				 //new_animal->previous_node=NULL;
				 animal_node_pointer[current_akey]->previous_node=NULL;
				 //animal_node_pointer[current_akey] = new_animal;
				// if(new_animal->akey ==2454078)
               // {
				
               // printf("previous is %lld next is %lld", new_animal->previous_node, new_animal->next_node);
          // }
				 }
				 
               }
           }
    /*Adding movement to event node is done*/
        //   printf("adding movement done");
        //   system("pause");
           
           

/* 2.ADD BIRTH EVENT AND ADD BIRTH ANIMALS To THE ANIMAL POINTER*/
    /* IGNORE //@@@
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
                
               
                //ADD THE NEW NODE TO THE ARRAY 
                 add_event_node(event_day, current_day, new_event) ;
                 
                
             //THEN CREATE ANIMAL NODE
                //struct animal_node *new_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
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
          @@@ */
/* 3.=========================ADD TESTING EVENTS================================================================*/
/* ONLY IF disease is TB*/
		if(disease==0)
		{
			    
           for (i=0; i < num_tests; i++)
          { 
                int current_day = (int)TestData[i][0] ;
                if(current_day<sim_days)
                {
				struct event_node *new_event;
                new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
                
                new_event -> event_type = 5;
                new_event -> src_pro_id = (long long)TestData[i][1] ; //this is testing schedule
				new_event -> des_pro_id = -100; // destination farm
                new_event -> src_testarea = -100 ;
				new_event -> next_node = NULL;   
         
                
               
                //ADD THE NEW NODE TO THE ARRAY 
                 add_event_node(event_day, current_day, new_event) ;
                }
           }
/*test adding done*/
	   } 
/*============================ADD TESTING EVENTS DONE=============================================================*/	   

/*++++++++++++++++++++++++++++++++++++++++++++++++++IN PROGRESS+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 4. ADD SAMPLING EVENTS*/
/*The simplest case is that sampling occurs at the same date.
But this might be too simplified for FMD if virus is already cleared from the host.
Set a regular interval for sampling?*/
/*Have to refresh memories on sampling in birth-death and coalescent models.
But the sampling should be random I believe. So for now assume random selection, which means sampling should be added
overtime and not as a priori.


Added 18 Sep 2017: Hall's paper recommends sampling should be uniform over time and location without a reference to effective population size or the number of infected individuals in that period/locations.
That means perhaps the best is to devide NZ into regions and then regularly collect samples out of each region with a given interval. 
Consider: if a particular region does not have any infected individuals, do I still collect a given number of total samples collecting more from other regions? or just do with whatever available?
Maybe the latter.

Sampling is broken down to two-fold:
(1) Control/detection strategy that actually collects samples
(2) Sample selection from those collected

Hall's paper is more about (2), but that ignores the fact that collected samples are actually already biased. This bias happens because sampling effort is higher in higher risk area or urban area.
Contact tracing also introduces bias in generalisability of samples. Important point here is that control/detection can itself change disease spread.

So, this study should focus on (1), but touch on (2) too.

Different control/detection strategies considered:
(a) uniform detection: ideal world? can detect randomly stratified to time and geographical locations.
(b) detection with different intensity: like bTB surveillance in New Zealand.
(c) detection with contact tracing: maybe uniform detection + contact tracing?

Sampling selection considered:
(A) uniform selection
(B) proportional to the number of samples in time and space
(C) selection based on community structure: say one sample from one community?

*/


/* everytime when update the infection pressure and so on, calculate the rate of detection.
This rate might be herd-level. And based on that rate, run the Gillespie to get the day for the next event.
So following scenarios go inside the function for markov event
*/




/*(a) Uniform detection*/
// so perhaps set a given probability for each herd to be detected and combine this into Gillespie algorithm?
// how to define this probability? Herd sensitivity based on tuberculin test? (that will vary depending on the number of infected animals in the herd...)
// maybe put some values and see how this affects the results? (sensitivity analysis on sensitivity)







/*(b) Detection with different intensity depending on regions*/
// difference can be due to frequency in testing, slaughterhouse capacity, or vets are more vigilant. set a different probability for each region.

/*(c) Uniform detection plus contact tracing*/
// in addition to a given probability, add backward and forward contact tracing
// how to define the detection mechanisms in traced herds? 
// Question: does contact tracing make population more homogeneous? or heterogeneous? or the same?


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		
           
/*================================ADD UNIT CHNAGE EVENTS=================================================================================*/
       int YEARS = 1; //@@@
           int current_day;
	    for(i=0;i<YEARS;i++)
           {
           	current_day = 365*i + 10;
           	//printf("current day is %d",current_day);
           	struct event_node *new_event; // I think just defining this once is enough (i.e. no need to include this in the loop)
                new_event = (struct event_node*)malloc(sizeof( struct event_node )); 
                
                new_event -> event_type = 8;
                new_event -> src_pro_id = -100 ; 
				new_event -> des_pro_id = -100; // destination farm
                new_event -> src_testarea = -100 ;
				new_event -> next_node = NULL;   
         
                
               
                /* ADD THE NEW NODE TO THE ARRAY */
                 add_event_node(event_day, current_day, new_event) ;
		   }
 /*================================ADD UNIT CHNAGE EVENTS DONE============================================================================*/
 
    
/*=========================DATA PREPARATION ENDS=============================================*/ 
//check if events were added properly
/* Following was used to check if events were properly added*/
	//printf("adding events done");
	//system("pause");	  
	  
 //for (i = 0 ; i < 12; i++)
   //   {
       visualize_list(event_day, 0);
     //     }
        //  printf("Before this you should have seen event lists!!");
//system("pause")	;	 
//  */


//visualize_animals(FarmProductionList,7186*3+1);
//system("pause");






/*=========================SIMULATION===================================================*/   
    

       
/* Below to check if non-Markov events in the event linked lists are working*/
/* TO DO*/
/*(1) NEEDS TO GIVE A DATE FOR NEXT NON-MARKOV EVENTS.
Checking each day if they have events and do these.
If reaches to the end of the day, then proceed to the next day when events occur.
Return the value for the date*/
double today_date = 1;
int next_non_markov_date;
double day_to_markov,updated_date;
long long move_akey, animal_akey;
double sum_inf_pressure;
struct event_node *current_event ;
//current_event = new_event = (struct event_node*)malloc(sizeof( struct event_node ));
struct event_node *adding_new_event ;

struct animal_node *current_animal ;
//int temp_days = 31;
int current_pro_id ;
int j = 0;
//int YEAR = 0;
//printf("Before 2055") ;
//visualize_animals(FarmProductionList,(2055*3+2)) ;
//printf("These are animals in 2056 before") ;
//visualize_animals(FarmProductionList,2056*3+2) ;
//system("pause");

//visualize_animals(FarmProductionList,1225*3+1);
//		system("pause");

//visualize_list(event_day,728);
//system("pause");
//	if(today_date==728)
//			{
//			for(i=1220;i<1226;i++)
//			{
  //          printf("id is %d", i);
	//		visualize_animals(FarmProductionList,i*3+1);
	//		system("pause");
	//	    }
	//	}
/*=======SIMULATION STARTS=================================================================================================================================================================*/	
while(today_date<sim_days)
{
	printf("this is day %f", today_date);
//	if(today_date==728||today_date==729)//727 ok 729 not ok, which means something happened on 728
//	{
//		for(i=1220;i<1226;i++)
//		{
//			printf("this is %d at day %d",i, today_date);
//		visualize_animals(FarmProductionList,i*3+1);
//		system("pause");	
//		}
//	}
	//{
	 //   printf("today is %d",today_date);
//		visualize_animals(FarmProductionList,1225*3+1);
//		system("pause");
//	}
	
	//printf("today is %d",today_date);
//	if(today_date==739)// 10, 375,.... between unit movement occurs
	//{printf("YEAR is now %d",YEAR);
	//	if(YEAR==2)
//		{
//	printf("yes");
//	visualize_animals(FarmProductionList,1225*3+1);
//	system("pause");
//		}
	//  move_production_type(FarmProductionStatus,FarmProductionList,num_total_farms);
	//  YEAR++;	
	//  printf("year is %d",YEAR);
	//}
  next_non_markov_date = ceil(today_date);
 // if(today_date>=720)
 // {
 // 	printf("that animal is in "animal_node_pointer[2454078]->current_pro_id
 // }
  while (event_day[next_non_markov_date] == NULL)
            {//loop2
		    next_non_markov_date++; // go to the next day
	    printf("go to next day") ;
            } // now get next non-markov date
  if(next_non_markov_date>=sim_days)
  {
  	break;
  }
  updated_date=update_markov_date(today_date,FarmData,FarmProductionStatus,FarmProductionList,num_farm_production,beta_a,beta_b,next_non_markov_date, event_day,
  detection_scenario, disease, mutation_model, max_S2toS3,seq_length, current_id_isolate, mu, nth_isolate, IsolateData, IsolateSeqData); //@need to change arguments
  printf("day is updated %f",updated_date) ;
// system("pause");
  if (updated_date==next_non_markov_date) // this means markov event did not happen
     {
	printf("yes updated date is next non markov date");
     //system("pause");
	 // loop 1
	  //this gets the first event on the day
     current_event = event_day[next_non_markov_date];
    // printf("current_event is %d",current_event->event_type);
     //printf("current event is %d",current_event->event_type);
     //system("pause");
     while(current_event!=NULL)
     {
     //	if(current_event->akey==1182045)
     //	{
     //		printf("this animal is from %lld to %lld with type %d",current_event->src_pro_id, current_event->des_pro_id,current_event->event_type);
	//	    printf("and current_pro_id is %lld and day %d",animal_node_pointer[1182045]->current_pro_id, today_date);
	//	    printf("previous animal is %lld and next animal is %lld",animal_node_pointer[1182045]->previous_node->akey,animal_node_pointer[1182045]->next_node->akey);
	//	 }
	 printf("current event is %d", current_event->event_type);
	   if (current_event-> event_type <= 4 ) // if movement or new birth or cull death
	      {
	        long long src_pro_id = current_event->src_pro_id;
	        long long src_farm_id = floor(src_pro_id/3) ;
	        int pre_move_det_status = FarmData[src_farm_id][6] ; //if detected 1, non detected 0
	        int testarea_num = current_event->src_testarea ;
	     //   if(Farm_ptr_inf_move[37761]==NULL)
	     //   {
	     //   	printf("Yes still NULL") ;
		//	}
		//	else{
		//		printf("No not any more NULL before") ;
		//		system("pause") ;
		//	}
			move_animal_unit(FarmProductionList,FarmData,FarmProductionStatus,current_event,animal_node_pointer,Se_occult,Se_detect, disease, Se_slaughter, today_date,
		    nth_isolate,nth_slaughter, inf_slaughtered, head_list_isolate, Farm_ptr_inf_move, IsolateData, IsolateSeqData
			); // function to move animals
		//	 if(Farm_ptr_inf_move[37761]==NULL)
	    //    {
	    //    	printf("Yes still NULL") ;
		//	}
		//	else{
		//		printf("No not any more NULL after") ;
		//		system("pause") ;
		//	}
			int after_move_det_status = FarmData[src_farm_id][6] ;
	    //    printf("movement done akey is %lld",current_event->akey);
	//	if(current_event->akey==1182045)
     //	{
     //		printf("this animal is from %lld to %lld with type %d",current_event->src_pro_id, current_event->des_pro_id,current_event->event_type);
	//	    printf("and current_pro_id is %lld and day %d",animal_node_pointer[1182045]->current_pro_id, today_date);
	//	    printf("and next animal is %lld",animal_node_pointer[1182045]->next_node->akey);
	//	 }
		  
	/*===========Here testing all animals in the herd plus contact tracing if infected animals are detected==========*/
	if(after_move_det_status == pre_move_det_status +1) // means this farm is detected in the movement process, either pre-movement or slaughter
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
         
     /* ADD THE NEW NODE TO THE ARRAY */
    add_event_node(event_day, today_date+5, new_event) ;
    }
	/*==========Movement section done==========*/
		
	      } // if movement done
	   else if (current_event-> event_type ==5 )//if this is testing
	      {
	   //  printf("testing");		
	 	   test_farms(FarmData,current_event,Se_occult,Se_detect) ; // get the testing schedule id
		 	 // if false-positive, they test only reacter animals?
		    //then no-FP + no-detection, no-FP+detection, FP+no-detection, FP+detection
	        //printf("test done") ;	
	        
	      } 
		    // if testing is happening, check all animals with test accuracy, then pick up test positive ones.
		    //Record information of positive animal (age, region of herd etc) and cull it.
		    //For positive herds, do follow up testing
	   else if (current_event-> event_type == 6||current_event-> event_type == 7) //if updating disease status
	      {
	      printf("updating disease status") ;
	      animal_akey = current_event-> akey ;
	        //current_animal = animal_node_pointer[animal_akey] ;
	        //current_pro_id = current_animal->current_pro_id;
	         current_pro_id = animal_node_pointer[animal_akey]->current_pro_id ;
	        int current_farm = (int)(floor(current_pro_id/3)) ;
	        printf("current farm is %d", current_farm) ;
		    if (current_event-> event_type ==6 ) //exposed to occult(TB) or exposed to infectiousness
		      {
		      	
		      //	printf("This is S2 to S3") ;
		      day_S3toS4 = rand()%(max_S3toS4)+1;
		     printf("day S3 to S4 is %d",day_S3toS4) ;
		      int day_to_add = day_S3toS4 + next_non_markov_date; //define when occult to detectable or infectious to recover happens
			  printf("day_to_add is %d",day_to_add) ;
				  	
		      animal_node_pointer[animal_akey]-> disease_status = 2 ;
		    //  printf("Num of status 2 is %f",FarmProductionStatus[current_pro_id][2] );
		      FarmProductionStatus[current_pro_id][2]--;
		    //  printf("Num of status 3 is %f",FarmProductionStatus[current_pro_id][3] );
		      FarmProductionStatus[current_pro_id][3]++;
		      FarmData[current_farm][7]++;
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
		    if (current_event-> event_type ==7) // occult to detectable (TB), infectious to immune (FMD)
		    //have to check how receovered are treated in coalescent/BD model
		      {
		      animal_node_pointer[animal_akey]-> disease_status = 3 ;
		      FarmProductionStatus[current_pro_id][3]--;
		      FarmProductionStatus[current_pro_id][4]++;
		      FarmData[current_farm][8]++ ;//increase detectable
		     // printf("detectable increased");
		     // system("pause");
		      FarmData[current_farm][7]-- ;//decrease occult
		      printf("FarmData updated") ;
		      //record isolate data
		      IsolateData[*nth_isolate][0] = animal_akey ;
		      IsolateData[*nth_isolate][1] = current_pro_id ;
		      IsolateData[*nth_isolate][2] = today_date;
		      IsolateData[*nth_isolate][3] = animal_node_pointer[animal_akey]->ptr_isolate->d_time;
	          IsolateData[*nth_isolate][4] = animal_node_pointer[animal_akey]->ptr_isolate->id;
		      IsolateData[*nth_isolate][5] = animal_node_pointer[animal_akey]->ptr_isolate->parent_id;
		      IsolateData[*nth_isolate][6] = animal_node_pointer[animal_akey]->ptr_isolate->num_mut;
			   IsolateData[*nth_isolate][7] = 0 ; // not isolated 
			   printf("IsolateData stored") ;
			  // IsolateSeqData[*nth_isolate][0] = current_animal->ptr_isolate->id ;
			   IsolateSeqData[*nth_isolate] = &(animal_node_pointer[animal_akey]->ptr_isolate->sequence);
			   *nth_isolate = *nth_isolate + 1 ;
			   printf("S3 to S4 done") ;
		      }
	       // printf("Updating TB status done") ;	
	      }//exposed->occult ot occult->detectable DONE
	      else if (current_event->event_type==8)
	      {
	      	printf("Changing unit!");
	      	move_production_type(FarmProductionStatus,FarmProductionList,num_total_farms); // @double check if I updated # of animals for each status?
	      //	printf("Changing unit done!");
	      //	system("pause");
		  }
	/*================Movement-ban lift==============================================================================*/
		  else if(current_event->event_type==10)
		  {
		  	 long long src_farm_id = floor(current_event->src_pro_id/3) ;
		  	 FarmData[src_farm_id][6] = 0 ; // farm detected status goes back to 0
		  }
	      
	      
	   // if(today_date==728)
     	//{
     	//	printf("this animal is from %lld to %lld with type %d",current_event->src_pro_id, current_event->des_pro_id,current_event->event_type);
		   // printf("and current_pro_id is %lld",animal_node_pointer[2454078]->current_pro_id);
		  // if(animal_node_pointer[1182045]->previous_node!=NULL&&animal_node_pointer[1182045]->previous_node->akey!=3546673)
		  // {
		  // 	if(animal_node_pointer[1182045]->previous_node!=NULL)
		   //	{
		   //	printf("previous animal is %lld",animal_node_pointer[1182045]->previous_node->akey);
		   	
		   //}
		    //if(FarmProductionList[1220*3+2]!=NULL)
		   //{
		   //	printf("1220 is pointing %lld",FarmProductionList[1220*3+2]->akey);
		  // }
		  // else
		  // {
		   //	printf("1220 is pointing null");
		  // }
		   	//printf("current akey is %lld", current_event->akey);
		   	//printf("And 1225 is pointing %lld",FarmProductionList[1225*3+1]->akey);
		   	//system("pause");
		   //}
		  // else if (animal_node_pointer[1182045]->previous_node==NULL)
		  // {
		   //	printf("1182045 is now first");
		   //	if(FarmProductionList[1220*3+2]!=NULL)
		  // {
		   //	printf("and 1220 is pointing %lld",FarmProductionList[1220*3+2]->akey);
		  // }
		  // }
		   
		// }  
	      
		struct event_node *previous_event;
		//previous_event =  (struct event_node*)malloc(sizeof( struct event_node )); // passing a pointer that points to a memory allocated dynamically
	   previous_event = current_event ; //rewire to the next event
	   current_event = current_event->next_node;
	   if (current_event!=NULL)
	   {
	 //  	printf("going to free event now") ;
	   free(previous_event);//this is not dynamic memory, but is it ok to free?
	   printf("next event is %d, %lld", current_event->event_type, current_event->akey);
       }
       else
       {
       	printf("next event is NULL") ;
	   }
	   event_day[next_non_markov_date] = current_event;
	   //printf("event day is pointing %lld",event_day[next_non_markov_date]->akey);
	  // j++;
	   //printf("this is %d th event",j);
	  // system("pause");
	 
	 /*---------------------------------------------------
	 Here add capturing the number of infected and do sampling.
	 Hold on: when randomly sampling, should I include all exposed ones or only infectious ones?
	 ---------------------------------------------------------------------------------------------*/ 
	  
     }//while loop for going to next events ends
	
	
     } // this is the end of loop 1    
     today_date = updated_date;
  }
  /*=============================SIMULATION ENDS==============================================================================================================================================================*/	
	
	
	



 


//printf("These are animals in 2055") ;
//visualize_animals(FarmProductionList,2055*3+2) ;
//printf("These are animals in 2056") ;
//visualize_animals(FarmProductionList,2056*3+2) ;	 
//for (i=0;i<num_total_farms;i++)
//{
//	if(FarmData[i][8]>0)
//	printf("this has detectable");
//  }  
	   
	   
count_farms_infected_detected(FarmData,OutPut,num_total_farms,iteration);
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
 //  free(FarmProductionStatus); // I don't have to free this memory either
   //free(FarmProductionList);
   
  // printf("Farm freed");
   
   /*Clear animal node*/
   // Do I need this as animal_node_pointer is just pointing to animal_node, which is deleted just above
   // Yes for now free only pointer to animal_node, not the double pointer
//  for(i=0; i < num_total_animals; i++)
 // {
   //free(animal_node_pointer[i]) ;
//   animal_node_pointer[i]->next_node = NULL;
//   animal_node_pointer[i]->previous_node = NULL;
//    }
   // free(animal_node_pointer) ;
   // printf("animal node fred");
   // free(event_day);
   // printf("event fred");
   
   export_IsolateData(OutIsolateDataFile, IsolateData, IsolateSeqData, num_max_isolate, num_isolate_vars) ;
   free(OutIsolateDataFile) ;
   
}//END OF EACH ITERATION

/*===========END OF ONE ITERATION==============================================================================================================================================*/


/*=======================EXPORT OUTPUT (END OF MAIN)========================================================================================================*/
write_OutPut(OutPutFile,OutPut,tot_iterations,num_OutPut) ;    
//@@@@ need to export isolate data; 
// NEEDS TO EXPORT THE RESULTS







} //END OF MAIN INT


/*======================EXPORT OUTPUT DONE=====================================================================================================================*/

	
	
	// NOW UPDATE THE EVENT TABLE SO THAT EVENTS THAT INVOLVE INFECTED FARMS ARE ONLY LISTED
	

	
	// then choose newly infected individuals randomly or following some rules
	// update their status - how to code choosing animals from arrayed structs?
	// then now have to trace the infection status of infected one to update their status 
	// when they become infectious from latent period.
	// Option1: define latent period first once they get infected and change status when it comes
	// Option2: Do Gllepise and set the date for the next event and when they occur choose event
	
	// Within-herd infection can be independently modelled unless movement occurs
	
	// Maybe the way is to look up the next global event which involves infected farm
	// then if the waiting time is after the global event, switch to the daily-basis event?
	// well I can't really think how to change the rate of events when global events occur.
	// So maybe better to update daily from beggining.
	// But because latent period is such long and transmission occurs so slow compared to
	// non-Markov events, maybe we just need to update at every non-Markov events.
	
	// Based on Andrew's note from workshop
	// Calculate waiting time and if this time is beyond the next non-Markov event
	// then do non-Markov and update time and get waiting time again.
	// else do Markov-event and update time and event.
	// Rate of event is the total rate of within-farm transmission and wildlife transmission
	// from all bTB infected and farms that can potentially infected from wildlife.
	// Just have to ask Andrew if using this approach is ok or not.
	
	
	
	// so do a bernouil trial using beta*S*I + beta_wildlife*S
	// Does each management group have different infection rate from wildlife?
	// Yes it should
	
	// repeat update everyday until some of the global events occur
	// Global events here include
	// on- and off-farm movement from bTB infected farm
	// culling and testing at infected farms
	// calving at the infected farms
	
/*3.2 Modelling management decision*/

	// decide when introduction and selling occur as a function of some farm characteristics
	// like high-turn over farm and trading farm
	// or auction date (which can be defined as priori because cattle move on the same date
	// so the precise date doesn't really matter
	// how about internet trading people? That would be much dependent on availability?
	// i.e. driven by the behaiour of source farm and traders
	
	// between management group occurs mandatory according to their age or months
	
	// need to know if move-in or move-out occurred on a given date.
	// update the herd size, number of infected and stuff
	// and go the next day
	
	// calving depends on calvng rate, which is a function of what?
	// does vertical transmission occurs if dam is btb infected and milk their calves?
	// culling depends on calving results
	


	
                
                
               
 
 
 
/*=====FUNCTIONS===============================================================================================================================*/


/* -------------------------------------------------------------------------- */
/* read_farm_data: READING AND PARSING CSV FARM LIST */
/* -------------------------------------------------------------------------- */
void read_farm_data(char FarmDataFile[], long long **FarmData,int num_farms)
{     
    /* OPEN INPUT FILE */
    FILE *Farms = fopen(FarmDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    long long line_num, farm_id, testarea, disease_status,test_schedule_id;
    double x_coord, y_coord;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_farms; line_num++)
      { 
         fscanf(Farms, "%lld,%lf,%lf,%lld,%lld, %lld", &farm_id, &x_coord, &y_coord, &testarea,&disease_status, &test_schedule_id);

         /* STORE VALUES IN FARM LIST */
             FarmData[line_num][0] = farm_id;
             FarmData[line_num][1] = x_coord;
             FarmData[line_num][2] = y_coord;
             FarmData[line_num][3] = testarea ;
             FarmData[line_num][4] = disease_status; 
             FarmData[line_num][5] = test_schedule_id;
             FarmData[line_num][6] = 0;//tb_detected or not
             FarmData[line_num][7] = 0;//sum of occult
             FarmData[line_num][8] = 0;//sum of detectable
             FarmData[line_num][9] = -1; //farm id for that will infect this farm first
             FarmData[line_num][10] = -1; // time that this farm gets infection from the above farm
             
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Farms);
   
} 
/* -------------------------------------------------------------------------- */
/* READ INITIAL ANIMAL DATA*/
void read_animal_data(char AnimalDataFile[], long long **AnimalData, int num_animals)
{
	
/* OPEN INPUT FILE */
    FILE *Animals = fopen(AnimalDataFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, farm_id, age_day,Type, Sex, Breed, pregnant_status, disease_status;
    long long akey;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_animals; line_num++)
      { 
         fscanf(Animals, "%lld,%d,%d, %d, %d,%d,%d,%d",&akey, &farm_id, &age_day, &Type, &Sex, &Breed,&pregnant_status,&disease_status );

         /* STORE VALUES IN AnimalData */
             AnimalData[line_num][0] = akey;
             AnimalData[line_num][1] = farm_id;
             AnimalData[line_num][2] = age_day;
             AnimalData[line_num][3] = Type ;
             AnimalData[line_num][4] = Sex; 
             AnimalData[line_num][5] = Breed;
             AnimalData[line_num][6] = pregnant_status;
             AnimalData[line_num][7] = disease_status; 
             //AnimalData[line_num][8] = AnimalData[line_num][1]*3+AnimalData[line_num][3];
              
             
      }            
   
   /* CLOSE INPUT FILE */
   fclose(Animals);
   printf("First animal is %lld",AnimalData[0][0]) ;
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
/*--------------------------------------------------------------------------------*/

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
			 printf("%d,%lld ", current_node1 -> event_type, current_node1 -> akey);
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
void move_animal_unit(struct animal_node *FarmProductionList[],long long **FarmData,double **FarmProductionStatus, struct event_node *current_event, 
struct animal_node **animal_node_pointer,double Se_occult, double Se_detect, int disease, double Se_slaughter, double date, int *nth_isolate, int *nth_slaughter, 
struct animal_node *inf_slaughtered[], struct isolate *head_list_isolate[], struct inf_move *Farm_ptr_inf_move[], double **IsolateData, char *IsolateSeqData[])
{
	int today_date = (int)date;
	int current_event_type = current_event->event_type;
	int src_pro_id = current_event->src_pro_id ;
//	printf("src_pro_id is %lld", src_pro_id) ;
	double random_value;
	//printf("src_pro_id is %d", src_pro_id);
	long long current_akey = current_event->akey ;
//	printf("akey is %lld", current_akey);
	int stop = 0;
	struct animal_node *moving_animal;
	
		//moving_animal = (struct animal_node*)malloc(sizeof( struct animal_node )); This is NOT right I think as this creates both a memory for the object as well as a pointer: we only need a memory for the pointer.
	moving_animal = animal_node_pointer[current_akey] ;
//	if(animal_node_pointer[current_akey]==NULL)
//	{
//		printf("this animal is null");
//		system("pause") ;
//	}
//	if(moving_animal->akey==0)
//	{
//		printf("akey is %lld", current_akey);
//		system("pause") ;
//	}
	
	int disease_status = moving_animal->disease_status ;
	//struct animal_node *added_animal;
	//added_animal = (struct animal_node*)malloc(sizeof( struct animal_node )); 
//	if(current_akey==2454078)
//	{
//		printf("previous node is %lld and next node is %lld",moving_animal->previous_node, moving_animal->next_node);
//	}
	
	
 /*========Movement to slaughter==============================================================================================================*/	
 if(current_event_type==4) //slaughter/death
 {
 	if(disease_status>=2&&disease==0) //if TB and the animal is occult(reacive) or detectable(infectious). No slaughter inspection for FMD.
 	{
 		double random_value = ((double)rand()/(double)RAND_MAX) ;
 		if(random_value<=Se_slaughter)
 		{
 		// change the herd status to infected
 		FarmData[src_pro_id][4]=1 ; //becomes infected (should be marked as infected already)
 		FarmData[src_pro_id][6]=1; // becomes detected 
 		// change the number of animals in the disease status group
 		FarmData[src_pro_id][disease_status+5] --; //minus 1
 		FarmProductionStatus[src_pro_id][disease_status+1] --; //minus 1
 		
/*================Remove animals from the linked lists=============================================================*/
 		      if (moving_animal -> previous_node != NULL)
	          { // if moving animal's previous node is conencted to other animal
	    //    
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
		        moving_animal->current_pro_id = -10 ; //means it's sent to slaughterhouse 
 		
 /*===============Removing the animal from the linked list done=================================================================================*/		
 
 /*===============Add the animal to infected-slaughtered linked list=============================================================================*/
 		if(inf_slaughtered[today_date]==NULL)
				{
				inf_slaughtered[today_date] = moving_animal; //
				}
		else
				{
				struct animal_node *next_slaughtered; 
				next_slaughtered = inf_slaughtered[today_date] ;
				inf_slaughtered[today_date] = moving_animal ;
				moving_animal->next_node = next_slaughtered ;
				moving_animal->previous_node = NULL ; // no need to record a previous animal anymore
				}
				
 		//collect sequence data
 		struct isolate *ptr_current_isolate; 
 		int nth = *nth_isolate ;
		head_list_isolate[nth] = moving_animal->ptr_isolate ; //handover the isolated data @@@ have to define nth_isolate and include it in the function argument
		//@@@@ why did I want this head_list_isolate? 
		ptr_current_isolate = moving_animal->ptr_isolate ;
		IsolateData[nth][0] = moving_animal->akey ;
		IsolateData[nth][1] = moving_animal->current_pro_id ;
		IsolateData[nth][2] = today_date;
		IsolateData[nth][3] = ptr_current_isolate->d_time;
		IsolateData[nth][4] = ptr_current_isolate->id;
		IsolateData[nth][5] = ptr_current_isolate->parent_id;
		IsolateData[nth][6] = ptr_current_isolate->num_mut;
		IsolateData[nth][7] = 1 ;
		
		//IsolateSeqData[nth][0] = moving_animal->ptr_isolate->id ;
		IsolateSeqData[nth] = &(moving_animal->ptr_isolate->sequence);
		//IsolateData[nth][8] = current_animal->ptr_isolate->sequence ;
		*nth_isolate = nth + 1;  // increment isolate counter by 1
		 	 
		 // testing all animals in the herd - this will be in the main function
 		// movement ban - this already done
 		// contact tracing - this will be in the main function
 		
		 }
	 }
	 else 
	 {
	 	*nth_slaughter = *nth_slaughter + 1 ;
	 }
 }
	
/*===========================Movement to slaughter END====================================================================================*/

/*===========================Movement except new calves born===============================================================================*/	
  else if (current_event_type != 3) // if not 3 it means it is not new born
  {
    //printf("Not new calve");  
	int src_farm_id = floor(src_pro_id/3) ;
	//printf("farm id is %d", src_farm_id);
/* Check if movement is allowed*/	
	if (current_event_type<=2 & FarmData[src_farm_id][6]==1) // if Disease status is detected here, cancel
	{
		stop == 1;
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
				FarmData[src_farm_id][6]=1 ; //farm becomes detected
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
		stop == 0;
	}


/*===================IF MOVEMENTS ALLOWED=================================================================================*/
if (stop == 0) // if the movement is still allowed
{


/* Decide if the moving animal is new or not*/
    //if (animal_pointer_list[current_akey]==NULL)
//	{
	//	printf("this is new animal");
/*If this animal is new*/
		
		//printf("current akey is %lld",current_akey);
		
		//added_animal->age_day = 0;
		//added_animal->akey = current_akey;
		//printf("b");
		//added_animal->type = current_event_type ;
		//printf("c");
		//added_animal->tb_status = 0 ;
		//printf("d");
		//added_animal->breed = 6 ;
		//printf("e"); 
		//added_animal->num_births = 0;
		//printf("f");
		//added_animal->sex = 1;
		//printf("g");
		//added_animal->pregnant_status = 0;
		//added_animal->current_pro_id = current_event->des_pro_id;
		//printf("h");
		//added_animal->next_node = NULL;
		//printf("i");
		//added_animal->previous_node = NULL;
		//printf("new animal node is made");
		//animal_pointer_list[current_akey] = added_animal;
		
		//printf("new animal is added");
	//}
/*If this animal is new end*/
	//else
	
/*If this animal is already existing*/
	//{//LOOPB
	//printf("this is not new animal");
	//if (FarmProductionList[src_pro_id] != NULL) // double check if there are animals
      // {//LOOPA
  	
  	/*===========Add this move to Farm_ptr_inf_move if the source farm is infected and if this move is not yet recorded*/
  	if(FarmData[src_farm_id][4]==1) // if this farm is infected
  	{
  	//	printf("This %lld farm is infected",src_farm_id) ;
  		long long des_farm_id = floor(current_event->des_pro_id/3);
  		if(FarmData[des_farm_id][9]==-1)
  		{
  			FarmData[des_farm_id][9] = src_farm_id; // record who infected this farm
		  }
  		if(Farm_ptr_inf_move[src_farm_id]==NULL)
  		{
  		//	printf("src_farm_id is NULL") ;
  		struct inf_move *new_node = (struct inf_move*)malloc(sizeof(struct inf_move));
		  new_node->date = today_date ;
		  new_node->des_farm_id = des_farm_id ;
		  new_node->next_node = NULL;
		  Farm_ptr_inf_move[src_farm_id] = new_node;
		  }
		else 
		{
	/*========NEEDS CHANGE IF INFECTED TO BE RECORDED============================*/		
		//	printf("src_farm_id %lld is not NULL", src_farm_id) ;
			
  // printf("The value of des is : %lld\n", Farm_ptr_inf_move[src_farm_id]->des_farm_id  );
  // system("pause") ;
		struct inf_move *next_node = (struct inf_move*)malloc(sizeof(struct inf_move));
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
		}
	/*========NEEDS CHANGE IF INFECTED TO BE RECORDED DONE=======================*/
	  } // recording infected movements done
  	       
	       
	   //  printf("this animal is %lld", moving_animal->akey) ;
	       if (moving_animal -> previous_node != NULL)
	          { // if moving animal's previous node is conencted to other animal
	    //      printf("A starts");
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
	   FarmProductionStatus[src_pro_id][0] = FarmProductionStatus[src_pro_id][0] - 1;
	   int current_status = moving_animal->disease_status ;
	   FarmProductionStatus[src_pro_id][current_status+1] = FarmProductionStatus[src_pro_id][current_status+1] - 1;	
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
	int current_status = moving_animal->disease_status ;
	FarmProductionStatus[des_pro_id][current_status+1]++;
	moving_animal->current_pro_id=des_pro_id;
	
   
} // if stop==0 ENDS

/*==========================IF MOVEMENTS ALLOWED ENDS===========================================================*/

} //if (current_event_type != 3) ENDS

else if (current_event_type == 3) // if this is new born just add
{
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
	moving_animal->current_pro_id=des_pro_id;
} // new born ENDS

}
    
    
/* -------------------------------------------------------------------------- */
/* TESTING HERDS */
/* -------------------------------------------------------------------------- */
void test_farms(long long **FarmData,struct event_node *current_event,double Se_occult, double Se_detect)
{
	int testschedule_id = current_event->src_pro_id;//NOTE: THIS IS NOT SRC ID, BUT TEST ID
	for (i=0; i<num_total_farms; i++)
	{
		if(FarmData[i][5]==testschedule_id) // if the farm has the same testing id
		{
			int num_occult = FarmData[i][7];
			int num_detectable = FarmData[i][8] ;
		  
		  if((num_occult+num_detectable)>0) // if this farm has more than 0 occult+detectable
		  {
		     double P_miss = (pow(1-Se_occult,num_occult))*(pow(1-Se_detect,num_detectable)) ;
			 if ( ((double)rand()/(double)RAND_MAX) > P_miss)
			  {
			  FarmData[i][6]=1 ; //farm becomes detected
				// now omit the function which animals to be detected
				//in future to add the detected animals to "detected list"
				
			  }
		  }	
		}
	}
}
/*----------------------------------------------------------------------------------*/
/* MOVING ANIMALS BETEWEEN PRODUCTION TYPE*/
/*----------------------------------------------------------------------------------*/
void move_production_type(double **FarmProductionStatus,struct animal_node *FarmProductionList[],int num_total_farms2)
{
struct animal_node *current_animal;
struct animal_node *next_animal;
struct animal_node *target_animal;


int current_status;

 for (i=0; i<num_total_farms;i++)
     {
    	//printf("%d",i);
     	while(FarmProductionList[i*3+1]!=NULL)//if there are heifer
     	{
     //	if(i>=1210 && i <= 1230&& FarmProductionList[i*3+1]->next_node!= NULL)
	//	    {
		    //	printf("doing %d", i);
		    //	printf("this animal is at %lld",animal)
		   // current_animal= FarmProductionList[i*3+1];
		   // if(current_animal->akey==2454078)
		    // {
		    // 	printf("Yes I am at %d",i) ;
			 //}
	//	    }
		//printf("not null %lld I is %d",FarmProductionList[i*3+1]->akey,i);
		//printf("next node is %lld",FarmProductionList[i*3+1]->next_node->akey);
		// system("pause");
	     //}
     //	current_animal = (struct animal_node*)malloc(sizeof( struct animal_node )); 
     //	next_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
     //	target_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
     	
		 current_animal = FarmProductionList[i*3+1]; //get the first animal
		 //if checking age, do here in future
	//	if(current_animal->akey==1182045)
	//	{
	//		printf("unit is %lld i is %d", current_animal->current_pro_id, i);
			//printf("this animal is now farm %lld and previous animal is %lld next is %lld",current_animal->current_pro_id,current_animal->previous_node->akey,current_animal->next_node->akey);
	//		if(current_animal->next_node!=NULL)
	//		{
	//			printf("next is %lld",current_animal->next_node->akey);
	//		}
		//	system("pause");
	//	}
		//this check indicates that
		current_animal->type++; 
		current_animal->current_pro_id++;
		 //adding the next animal to farm pointer
		 
		
		 if(current_animal->next_node!=NULL)
		 {
		 next_animal = current_animal->next_node ;
		 next_animal->previous_node=NULL;
		 FarmProductionList[i*3+1] = next_animal;
	     }
	     else
	     {
	     	FarmProductionList[i*3+1] = NULL;
		 }
		 current_animal->previous_node = NULL;
		 current_animal->next_node = NULL;
		 
		 //update number of animals in the unit
		 current_status = current_animal->disease_status;
		 FarmProductionStatus[i*3+1][current_status+1]--;
		 FarmProductionStatus[i*3+1][0]--;
		 
		 
		 //adding the current animal to the destination
		  // first animal in adult list
		 if(FarmProductionList[i*3+2]==NULL)
		 {
		 	FarmProductionList[i*3+2] = current_animal;
		 }
		 else
		 {
		 	target_animal = FarmProductionList[i*3+2];
		 	target_animal->previous_node = current_animal;
		 	FarmProductionList[i*3+2] = current_animal;
		 	current_animal->next_node = target_animal;
		 }
		 
     	 FarmProductionStatus[i*3+2][current_status+1]++;
		 FarmProductionStatus[i*3+2][0]++;
		 
	    }//while ends
  //   	printf("heifer ends");
     	//next calf
     	while(FarmProductionList[i*3]!=NULL)//if there are heifer
     	{
     		
     	//current_animal = (struct animal_node*)malloc(sizeof( struct animal_node )); 
     	//next_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
     	//target_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
     	
		 current_animal = FarmProductionList[i*3]; //get the first animal
		 //if checking age, do here in future
		current_animal->type++;
		current_animal->current_pro_id++; 
	//	if(current_animal->akey==2454078)
	//	{
	//		printf("calf unit is %lld i is %d", current_animal->current_pro_id, i);
			//printf("this animal is now farm %lld and previous animal is %lld next is %lld",current_animal->current_pro_id,current_animal->previous_node->akey,current_animal->next_node->akey);
	//		if(current_animal->next_node!=NULL)
	//		{
	//			printf("next is %lld",current_animal->next_node->akey);
	//		}
			//system("pause");
	//	}
		 //adding the next animal to farm pointer
		 
		 if(current_animal->next_node!=NULL)
		 {
		 	next_animal = current_animal->next_node ;
		 next_animal->previous_node=NULL;
		 FarmProductionList[i*3] = next_animal;
	     }
	     else
	     {
	     	FarmProductionList[i*3] = NULL ;
		 }
		 current_animal->previous_node = NULL;
		 current_animal->next_node = NULL;
		 
		 //update number of animals in the unit
		 current_status = current_animal->disease_status;
		 FarmProductionStatus[i*3][current_status+1]--;
		 FarmProductionStatus[i*3][0]--;
		 
		 
		 //adding the current animal to the destination
		  // first animal in adult list
		 if(FarmProductionList[i*3+1]==NULL)
		 {
		 	FarmProductionList[i*3+1] = current_animal;
		 }
		 else
		 {
		 	target_animal = FarmProductionList[i*3+1];
		 	target_animal->previous_node = current_animal;
		 	FarmProductionList[i*3+1] = current_animal;
		 	current_animal->next_node = target_animal;
		 }
		 
     	 FarmProductionStatus[i*3+1][current_status+1]++;
		 FarmProductionStatus[i*3+1][0]++;
		 
	    }//calf ends
//	    printf("calf ends %d",i);
	 }  	
	
//printf("now free memory");
  //       free(current_animal);
    //     free(target_animal);
    //     free(next_animal);
//printf("ends");	
}




/*----------------------------------------------------------------------------------*/
/* CALCULATE NEXT MARKOV EVENT DATE*/
/*----------------------------------------------------------------------------------*/
double update_markov_date(double today_date, long long **FarmData, double **FarmProductionStatus,struct animal_node *FarmProductionList[],int num_farm_production,double beta_a,double beta_b, 
int next_non_markov_date, struct event_node *event_day[], int detection_scenario, int disease, int mutation_model, int max_S2toS3, int seq_length, 
int *current_id_isolate, double mu, int *nth_isolate, double **IsolateData, char *IsolateSeqData[])
{

double inf_pressure, mutation_pressure,inf_pressure_wild,sum_pressure,cumu_pressure, detection_pressure;
sum_pressure = 0;
cumu_pressure = 0;
detection_pressure = 0;
double day_to_markov;
int farm_id ;
int k = 0;
int testareanum ;

double tb_detection_p;
tb_detection_p=0.9; //@@@@do I put higher sensitivity for this?

 for(i=0; i<num_farm_production; i++)
	{//for loop A
	
	farm_id = (int)floor(i/3) ;
	testareanum = FarmData[farm_id][3] ;
	
	/*=================Calculate the total rate of event in Gillespie========================================================*/
	if ((testareanum==0||FarmProductionStatus[i][3]>0||FarmProductionStatus[i][4]>0||FarmProductionStatus[i][2]>0)&&FarmProductionStatus[i][0]>0)
	{    //either: (farm is in movement control || there is occult animal || there is detectable animal) AND there is susceptible animal AND there is animal in the herd
		//if(FarmProductionStatus[i][0]<0)
		//{
		//	printf("N is <0 something wrong!");
		//}
		//printf("yes this farm is %d",i);
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
			if(mutation_model==1) 
			{
				mutation_pressure = mu*(FarmProductionStatus[i][2]+FarmProductionStatus[i][3]+FarmProductionStatus[i][4]) ; //mutation occurs in exposed, occult, and detectable equally
				FarmProductionStatus[i][7] = mutation_pressure;
			} //@ set mu and mutation_pressure parameter
	    }
	    else if (disease==1)//FMD
	    {
	    inf_pressure = beta_a*FarmProductionStatus[i][1]*(FarmProductionStatus[i][3])/FarmProductionStatus[i][0]  ; //bSI/N	
	    //if(FarmProductionStatus[i][3]>0)
	    //{
	    //printf("inf pressure is %lf",inf_pressure) ;
	   // system("pause") ;	
		//}
		
	    if(mutation_model==1) //@ add mutation_model parameter
			{
				mutation_pressure = mu*(FarmProductionStatus[i][2]+FarmProductionStatus[i][3]); //mutation occurs in exposed and infectious equally
				if(mutation_pressure>0)
				{
				printf("mut pressure is %lf",mutation_pressure) ;
	   			system("pause") ;
				}
				
				FarmProductionStatus[i][7] = mutation_pressure;
			}
		}
		FarmProductionStatus[i][5] = inf_pressure; //inf_pressure, stored in column 5, indicates the total infection pressure in the herd
	/*========Calculate detection pressure and add it to the total rate=============================================*/
		if(mutation_model==1) // this needs to be changed in future if detection pressure is considered when it's not mutation model
		{
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
			    if(FarmData[farm_id][6]==0&&FarmProductionStatus[i][3]>0) // this farm should not be already detected
				{
				detection_pressure = 0.3/365;
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
				
			} // not store detection pressure but can easily calculated by [6] - [7] - [5]: CORRECTION - THIS IS NOT TRUE. [6] is a cumulative pressure over herds up to the current herd, so this formula won't give detection_pressure.
			// Should I store detection_pressure or calculate it again given that it sounds useless to store such a simple number in the data?
			// for now calculate it again but can be changed later
			else
			{
					detection_pressure = 0 ;
			}
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
	//printf("inf pressure is %f",sum_pressure);
	 double random_value = ((double)rand()/(double)RAND_MAX) ;
	 printf("random value is %f",random_value);
     day_to_markov =(-log(random_value))/sum_pressure ; // Waiting time
     printf("day_to_markov is %f",day_to_markov);
    // system("pause") ;
//ASSESS IF NEXT MARKOV TIME IS BEFORE NEXT NON MARKOV TIME     
	if (next_non_markov_date>day_to_markov+today_date)
	{ // if markov comes first, choose markov events 
	   k++;
	   //printf("This is %d th markov",k);
	   today_date = day_to_markov+today_date ;
	   //printf("ceil days is %d",today_date);
	   //system("pause");
	   double random_value = (double)rand()/(double)(RAND_MAX)*sum_pressure;
	   
/*=======Define the farm on which an event will occur====================*/
	   
	   
	   double current_value = 0;
	   int pro_id = 0;
	   
	   while(random_value>current_value) //random_value defines which farm event will occur at
	   {
	   	
	   	if(mutation_model==1)
	   		{
	   		if(FarmProductionStatus[pro_id][6]>0) //either infection or mutation pressure has to be >0 for the farm to have events
	   		
			   {
			   current_value = FarmProductionStatus[pro_id][6] ;
	   		printf("current value is %f",current_value) ;
	   		printf("random_value is %f",random_value);
	   		printf("sum_pressure is %f",sum_pressure);
	   	        }
		    }
		else if(mutation_model==0)
			{
			if(FarmProductionStatus[pro_id][5]>0)
			current_value = FarmProductionStatus[pro_id][6] ;
			}
		if(current_value>=random_value) // means the farm that has sum_pressure closest to but not bigger than random_value will be selected as the farm where an event occur
	       {
	       	printf("farm determined") ;
			   break;
	    		   }
		 else
		 {
		 	pro_id++;	
		 printf("pro id is %d", pro_id);
		 }    
		   //(disease==0&&(FarmData[farm_id][3]==0||FarmProductionStatus[pro_id][3]>0||FarmProductionStatus[pro_id][4]>0)&&(FarmProductionStatus[pro_id][0]!=0)) //the last one, it used to be FarmProductionStatus[pro_id][0] means if this is equal to 0, which should not happen
		   //FarmData[farm_id][3] is testarea, so checking if this is MCA
		   //next one evaluates if occult exists, next one for detectable exists
		   //ok, now think about FMD
	     
	     //  	printf("current_value is %f",current_value);
	      // 	printf("random value is %f",random_value);
	       
	   
	   }//pro_id is the farm to choose for the event
	   
	   farm_id = (int)floor(pro_id/3) ; // where farm_id was used??
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
		printf("mut model") ;
	// here to calculate the detection pressure again
		if(detection_scenario==0&&disease==1)
		{
		if(FarmData[farm_id][6]==0&&FarmProductionStatus[pro_id][3]>0)
			{
			detection_pressure = 0.3/365 ;
			}
			else
			{
			detection_pressure = 0;
			}
		printf("detection_pressure is %f",detection_pressure) ;
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
	double random_value = (double)rand()/(double)RAND_MAX*(FarmProductionStatus[pro_id][7]+FarmProductionStatus[pro_id][5]+detection_pressure); // get random value between 0 and sum of inf/mutation/detection pressure
    }
    else
    {
    double random_value = (double)rand()/(double)RAND_MAX*(FarmProductionStatus[pro_id][5]);
	}
/*=====================CHOOSE WHICH EVENTS?==============================================================*/	
printf("random value 2 is %f", random_value) ;
/*===================================WITHIN HERD TRANSMISSION PART STARTS===========================================================================================*/
	if(random_value<=FarmProductionStatus[pro_id][5])
	{
		printf("This is transmission") ;
	//	system("pause") ;
	//now get the first animal in the production unit
	   struct animal_node *current_animal;
	   //current_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
	   current_animal = FarmProductionList[pro_id];
	   int total_s = FarmProductionStatus[pro_id][1]; // number of susceptible	
	   int random_n = rand()%total_s+1;
	   int counter;
	   counter = 0;
	   while(current_animal!=NULL)
	    {//choosing which animal to infect begins
	    if(current_animal->disease_status==0) // if susceptible
	    {
	    counter++;	
		}
		if(counter==random_n)
			{ // add passing sequence here
			current_animal->disease_status=1;
			FarmProductionStatus[pro_id][1]--; //susceptible - 1
			FarmProductionStatus[pro_id][2]++; //Exposed + 1
			
			break;
	    	}
	    else
	     {
	     current_animal= current_animal->next_node ;		
		 }
	    
	    } // picked up an individual that will be infected
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
		random_n = rand()%num_infectious + 1; // a pseudo-random number between 1 and num_infectious
		counter = 0 ;
		struct animal_node *infecting_animal ;
		infecting_animal = FarmProductionList[pro_id];
	
	while(infecting_animal!=NULL)
	{
		if((infecting_animal->disease_status==2&&disease==1)||(infecting_animal->disease_status==3&&disease==0))// FMD and this animal is infectious
		{
			counter++;
		}
		if(counter==random_n)
		{
			//copy seq
			// ptr_isolate is pointing NULL, no memory allocated.
			struct isolate *new_isolate;
			new_isolate = (struct isolate*)malloc(sizeof( struct isolate ));
			strcpy(new_isolate->sequence, infecting_animal->ptr_isolate->sequence) ;// copy sequence
			new_isolate->d_time = today_date;
			new_isolate->parent_id = infecting_animal->ptr_isolate->id ;
			*current_id_isolate = *current_id_isolate + 1 ; // increment isolate id before assigning: this should be done for other occasions to be consistent
			new_isolate->id = *current_id_isolate; 
			current_animal->ptr_isolate = new_isolate ; // assign the pointer
		//	printf("New infection is %lld",infecting_animal->akey) ;
		//	system("pause") ;
			break;	
		}
		else
			{
			infecting_animal= infecting_animal->next_node ;	
			}
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
			  
			  add_event_node(event_day, day_to_add, adding_new_event) ; //am I using function here???? Inside function????
		      }
	
       	}
	      
	   
	/*===================================TRANSMISSION PART ENDS==============================================================================================*/


	else if(random_value>FarmProductionStatus[pro_id][5]&&random_value<=(FarmProductionStatus[pro_id][5]+FarmProductionStatus[pro_id][7])) // if mutation is selected
	/*==============MUTATION PART STARTS=================================================================================================================*/
		{
			printf("This is mutation") ;
		//	system("pause") ;
			//get random value n which is euql to the number of suscptible or that are subjcted to mutations
			//get nth animals that are eligible
		int num_tot_mutatable = (int)(FarmProductionStatus[pro_id][7]/mu) ; //total number of animals that are subjected to mutations
		int random_n = rand()%num_tot_mutatable+1; //get random number between 1 and num_tot_mutatable
		printf("A") ;
		 struct animal_node *current_animal;
		 printf("B") ;
		 current_animal = FarmProductionList[pro_id];
		 printf("C") ;
		 if(current_animal==NULL)
		 {
		 	printf("Something wrong! Animals should exist here but not!");
		 }
		 else
		 {
		 	printf("current animal is not null") ;
		 	int counter;
		 	counter = 0;
		 	while(counter<random_n)
		 	{
		 	if((disease==0&&current_animal->disease_status!=0)||(disease==1&&(current_animal->disease_status==1||current_animal->disease_status==2))) //TB disease==0, FMD 1. if TB then if animals are either exposed/
			 counter++; // increment
			 
			 if(counter==random_n)
			 {
			 	struct isolate *mutating_isolate ;
			 	mutating_isolate = current_animal->ptr_isolate; //now have to choose which base to change
			 	/*============================================================================================
			 	MUTATION MODEL HERE*/
			 	int r_position = rand()%seq_length ; //random position at nucleotide seq (random number generator needs improvement)
				char next_base; // make sure it is just char, not char[], which becomes character string type
				char current_base ; //get the current base
    			current_base = mutating_isolate->sequence[r_position] ;
    			printf("current base is %c", current_base) ;
    			if(current_base=='A') // make sure the character is bracketted with '', not " "
						{
							printf("Matched!");
							next_base= (char)"TGC"[rand()%3];
	
						}
				else if(current_base=='T')
						{
						printf("Matched!");
						next_base= (char)"AGC"[rand()%3];	
						}
				else if(current_base=='G')
						{
						printf("Matched!");
						next_base= (char)"ATC"[rand()%3];
						}
				else if(current_base=='C')
						{
						printf("Matched!");
						next_base= (char)"ATG"[rand()%3];	
						}
				else
						{
						printf("Error! This is not eligible base!") ;
						}
				printf("Next base is %c", next_base);
				mutating_isolate->sequence[r_position] = next_base;
				printf("Now base is %c", mutating_isolate->sequence[r_position]);
				mutating_isolate->num_mut++; //mutation count+
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
		 
		} // END OF MUTATION
	/*===================================MUTATION PART ENDS================================================================================================*/
	
	/*===================================DETECTION PART STARTS=======================================================================*/
	// if FMD assume all infected get identified and they have movement bans? (maybe not remove as many developing countries won't do so)
	// if no contact tracing then change the FarmProductionStatus[1]-[4]
	// and identify infected one through AnimalData dataframe, collect a random sample
	// how long for movement bans? In Argentina, by Perez (2004), movement restrictions lasted for more than 30 days, maybe do this?
	// after movement ban, all animals get completly disease clear status
	// Contact tracing: any herds that had contacts with the detected herds will have a movement ban: how about the secondary contact? Too much?
	
	// if bTB, within an infected herd, all or part of infected animals will be detected. 
	// two scenarios: one is complete removal in a defined time say 12 months after the movement ban.
	// The other is a realistic case where DSe is treated as inperfect and allows an escape of infected animals.
	else //else detection
	{
		printf("This is detection") ;
		system("pause") ;
		struct animal_node *current_animal;
		 current_animal = FarmProductionList[pro_id];
		 struct animal_node *previous_animal;
		 struct animal_node *next_animal ;
		 struct animal_node *next_slaughtered ;
		 
		 //struct isolate *current_isolate; 
		 if(current_animal==NULL)
		 {
		 	printf("Something wrong! Animals should exist here but not!");
		 }
		 if(FarmData[farm_id][6]==1)
		 {
		 	printf("This farm is already detected and should not be re-detected again!");
		 }
		 else //else1
		 {
		 	while(current_animal!=NULL) // until it reaches till the end
		 	{
		 	if(disease==1&&(current_animal->disease_status==2||current_animal->disease_status==1)) //if FMD infectious: for now if TB no passive detection occurs, if needed, add this: "(disease==0&&current_animal->disease_status>=2)||(" before hand and ")" in the end
			  if(detection_scenario==0)
			  {
			  	//random_value = ((double)rand()/(double)RAND_MAX) ; //get random value for determining if this animal would be detected or not
			  	//if(disease==0&&random_value>tb_detection_p) //this is the case where animals are not detected
			  	//	{
			  	 //     	current_animal = current_animal->next_node ; //next animal
			  	 //   }
			  	   	// first obtain seq data
			  	   		head_list_isolate[*nth_isolate] = current_animal->ptr_isolate ; //handover the isolated data @@@ have to define nth_isolate and include it in the function argument
			  			IsolateData[*nth_isolate][0] = current_animal->akey ;
			  			IsolateData[*nth_isolate][1] = current_animal->current_pro_id ;
			  			IsolateData[*nth_isolate][2] = today_date;
			  			IsolateData[*nth_isolate][3] = current_animal->ptr_isolate->d_time;
			  			IsolateData[*nth_isolate][4] = current_animal->ptr_isolate->id;
			  			IsolateData[*nth_isolate][5] = current_animal->ptr_isolate->parent_id;
			  			IsolateData[*nth_isolate][6] = current_animal->ptr_isolate->num_mut;
			  			IsolateData[*nth_isolate][7] = 1 ; //detected
			  			
					//	IsolateSeqData[*nth_isolate][0] = current_animal->ptr_isolate->id ;
					//strcpy(new_isolate->sequence, infecting_animal->ptr_isolate->sequence) ;/
			            IsolateSeqData[*nth_isolate] = &(current_animal->ptr_isolate->sequence);
						int i= 0 ;  
			            while(IsolateSeqData[*nth_isolate][i]!=NULL)
        				{
        					printf("%c",IsolateSeqData[*nth_isolate][i]) ;
						i++ ;
						}
						system("pause") ;
						  //IsolateData[*nth_isolate][8] = current_animal->ptr_isolate->sequence ;
			  			*nth_isolate = *nth_isolate + 1 ; // increment isolate counter by 1
			  			
			  			
			  	// Change the disease status and number recorded in FarmProductionStatus
			  			current_animal->disease_status= 3; // becomes immune - then have to cancel the event for them to beome immune or just skip it by saying only if no immune , immune occurs			  	
					if(current_animal->disease_status==1)
					{
						FarmProductionStatus[pro_id][2] -- ;
						FarmProductionStatus[pro_id][4]++ ;
						}	
			  		else if(current_animal->disease_status==2)
			  		{
			  			FarmProductionStatus[pro_id][3]--; //FMD fectious minus 1
			  			FarmProductionStatus[pro_id][4]++ ; //One increase in the recovered
					  }
			  			
			  			current_animal = current_animal->next_node ; //next animal
			  	}
			  } //while loop done
			  FarmData[farm_id][6] = 1 ; //becomes detected
			  	// and add a movement-ban-lift a month later
			  	struct event_node *adding_new_event;
			  			adding_new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		      			adding_new_event->akey=-1; // akey does not matter
		      			adding_new_event->src_pro_id = pro_id; 
		      			adding_new_event->des_pro_id = -100;
		      			adding_new_event->src_testarea=-100;
			  			adding_new_event->event_type=10;//occult to detectable happens
			  
			 			 add_event_node(event_day, floor(today_date)+30, adding_new_event) ;
		} // else1 done
			  		
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
void count_farms_infected_detected(long long **FarmData, double **OutPut, long num_total_farms2, int sim)
{
long i;
OutPut[sim][0]=0;
OutPut[sim][1]=0;
for (i=0;i<num_total_farms;i++)
{
	if(FarmData[i][7]>0||FarmData[i][8]>0)
	{
	 OutPut[sim][0]++; //infected	
	 //printf("one more infected");
	}
	if(FarmData[i][6]==1)
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
int export_IsolateData(char* OutIsolateDataFile, double **IsolateData, char *IsolateSeqData[], int num_max_isolate, int num_isolate_vars)
{
	FILE *Out = fopen(OutIsolateDataFile, "w") ;
	int line_num, col_num ;
	int i = 0;
	for (line_num = 0 ; line_num < num_max_isolate; line_num ++)
	{
		for (col_num = 0;col_num <num_isolate_vars ; col_num++ )
		
		{
	    fprintf(Out,"%lld,",IsolateData[line_num][col_num]);
        }
        while(IsolateSeqData[line_num][i]!=NULL)
        {
        	fprintf(Out,"%c",IsolateSeqData[line_num][i]) ;
			i++ ;
		}
        
        fprintf(Out,"\n");
    }
    fclose(Out) ;
    return 0 ;
	
}
