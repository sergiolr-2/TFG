#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
#include<time.h>
#include<sys/timeb.h>

#include <mpi.h>

struct row_group{// Every group have some value and belong to some sets
	int value;
	//int set_num[15];
};

struct row{ // A row is made of some 4bit group
	//int row_col;//4bit makes one column
	struct row_group *rg;//
};

struct cluster{ // A bicluster is made of some rows and some col
	int **matrix;
	int	rsize,csize;// row and col info
	struct row *rw;
	int row_col;
};

struct bicluster_component{ // Each member of bicluster is of this type....
	int rank;
	int discarded; // 0 or 1. Used for "lazy deleting"
	int row_pair1; // Tuple of rows that generated this cluster
	int row_pair2;
	int row_total;
	int col_total;
	int *row_number;
	int *column_number;
};

struct bicluster{   // bicluster that we can found....
	struct bicluster_component *bcc;
	int total_bicluster;
};

struct value{
	int *set;
	int total_set;// total number of sets
};

int binaryTodecimal(int r, int j);
void initializeSet();
void RUBIC();
int templateFound(int *);
int templateFound_inFile(int *,int); //,char *);
int search(int,int, int*);
void decimalTobinary(struct bicluster_component *b);
int binarySearch(int x,int t);
void decimalTobinary_fromTemplate(int *);
void patternConversionDecimalToBinary(int *t,struct bicluster *pattern);
void New_patternConversionDecimalToBinary(int *t,struct bicluster *pattern);
void printBicluster();
int New_templateFound(char *);
int mncRequirement(int *);
char* itxoa(int value, char* buffer, int base);

struct cluster clst;
struct bicluster bc;
struct value v[16];
int templ_num=0,**templ; // to templ means rho
FILE *fp2,*fp4,*fp5;
int mnr,mnc,num_bicluster=0;
int nprocs, rank; // MPI related variables


//----------------------------
// Makes process 0 receive all the clusters found
void gather_resuls(){
	// All processes flatten their clusters into one buffer.
	int local_bytes = 0;

	for (int i = 0; i < bc.total_bicluster; i++){
		local_bytes += sizeof(int) * (6 + bc.bcc[i].row_total + bc.bcc[i].col_total);
	}

	unsigned char* local_buffer = (unsigned char*) malloc(local_bytes);
	unsigned char* ptr = local_buffer;

	for (int i = 0; i < bc.total_bicluster; i++){
		// Single integers
		memcpy(ptr, &bc.bcc[i].rank, sizeof(int));
		ptr += sizeof(int);
		memcpy(ptr, &bc.bcc[i].discarded, sizeof(int));
		ptr += sizeof(int);
		memcpy(ptr, &bc.bcc[i].row_pair1, sizeof(int));
		ptr += sizeof(int);
		memcpy(ptr, &bc.bcc[i].row_pair2, sizeof(int));
		ptr += sizeof(int);
		memcpy(ptr, &bc.bcc[i].row_total, sizeof(int));
		ptr += sizeof(int);
		memcpy(ptr, &bc.bcc[i].col_total, sizeof(int));
		ptr += sizeof(int);
		// Arrays
		memcpy(ptr, bc.bcc[i].row_number, bc.bcc[i].row_total * sizeof(int));
		ptr += bc.bcc[i].row_total * sizeof(int);
		memcpy(ptr, bc.bcc[i].column_number, bc.bcc[i].col_total * sizeof(int));
		ptr += bc.bcc[i].col_total * sizeof(int);
	}

	// Process 0 obtains how many clusters each process will send and calculates the offsets betweeen each processes' buffer
	int* recv_bytes;
	int* offsets;
	unsigned char* global_buffer;
	int total_bytes = 0;

	if (rank == 0)
		recv_bytes = (int*) malloc(nprocs * sizeof(int));

	MPI_Gather(&local_bytes, 1, MPI_INT, recv_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0){
		offsets = (int*) malloc(nprocs * sizeof(int));
		for (int i = 0; i < nprocs; i++){
			offsets[i] = total_bytes;
			total_bytes += recv_bytes[i];
		}
		
		global_buffer = (unsigned char*) malloc(total_bytes);
		ptr = global_buffer;
	}

	// All processes send their buffer to process 0
	MPI_Gatherv(local_buffer, local_bytes, MPI_BYTE, global_buffer, recv_bytes, offsets, MPI_BYTE, 0, MPI_COMM_WORLD);

	// Get total number of biclusters (structs) received
	if (rank == 0)
		MPI_Reduce(MPI_IN_PLACE, &num_bicluster, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	else 
		MPI_Reduce(&num_bicluster, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	// Process 0 frees its previous bc structure and replaces it by the deserialized global buffer
	if (rank == 0){
		for (int i = 0; i < bc.total_bicluster; i++){
			free(bc.bcc[i].row_number);
			free(bc.bcc[i].column_number);
		}
		free(bc.bcc);

		bc.total_bicluster = num_bicluster;
		bc.bcc = (struct bicluster_component*) calloc(bc.total_bicluster, sizeof(struct bicluster_component));

		for (int i = 0; i < bc.total_bicluster; i++){
			//Single integers
			memcpy(&bc.bcc[i].rank, ptr, sizeof(int));
			ptr += sizeof(int);
			memcpy(&bc.bcc[i].discarded, ptr, sizeof(int));
			ptr += sizeof(int);
			memcpy(&bc.bcc[i].row_pair1, ptr, sizeof(int));
			ptr += sizeof(int);
			memcpy(&bc.bcc[i].row_pair2, ptr, sizeof(int));
			ptr += sizeof(int);
			memcpy(&bc.bcc[i].row_total, ptr, sizeof(int));
			ptr += sizeof(int);
			memcpy(&bc.bcc[i].col_total, ptr, sizeof(int));
			ptr += sizeof(int);

			//Arrays
			bc.bcc[i].row_number = (int*) calloc(bc.bcc[i].row_total, sizeof(int));
			bc.bcc[i].column_number = (int*) calloc(bc.bcc[i].col_total, sizeof(int));

			memcpy(bc.bcc[i].row_number, ptr, bc.bcc[i].row_total * sizeof(int));
			ptr += bc.bcc[i].row_total * sizeof(int);
			memcpy(bc.bcc[i].column_number, ptr, bc.bcc[i].col_total * sizeof(int));
			ptr += bc.bcc[i].col_total * sizeof(int);
		}
		// puts("Final result:");
		// for (int i = 0; i< bc.total_bicluster; i++){
		// 	printf("%d - rank=%d - discarded=%d - row_pair1=%d - row_pair2=%d - row_total=%d - col_total=%d\n",i, bc.bcc[i].rank, bc.bcc[i].discarded, bc.bcc[i].row_pair1, bc.bcc[i].row_pair2,bc.bcc[i].row_total, bc.bcc[i].col_total);
		// 	printf("r :");
		// 	for (int j = 0; j < bc.bcc[i].row_total; j++) printf("%d ", bc.bcc[i].row_number[j]);
		// 	printf("\nc :");
		// 	for (int j = 0; j < bc.bcc[i].col_total; j++) printf("%d ", bc.bcc[i].column_number[j]);
		// 	printf("\n");
		// }
	}

	// Free memory
	free(local_buffer);
	if (rank == 0){
		free(global_buffer);
		free(recv_bytes);
		free(offsets);
	}

}

//----------------------------
// Auxiliariy function for qsort to sort the bc struct according to the order in which each bicluster would appear in the original code (sequential)
int sequential_order(const void *a, const void *b){
	const struct bicluster_component* c1 = (const struct bicluster_component*) a;
	const struct bicluster_component* c2 = (const struct bicluster_component*) b;

	if (c1->row_pair1 != c2->row_pair1)
		return c1->row_pair1 - c2->row_pair1;
	else 
		return c1->row_pair2 - c2->row_pair2;
}

//----------------------------
// Algorithm that makes process 0 discard those clusters that wouldn't appear in the original code (sequential)
void discard_clusters(){
	if (nprocs == 1) return;

	int i, j, k, t, curr_div, last_div, tsize, found_component;
	// bc.bcc[i] represents the template that would enter templateFoundInFile
	// bc.bcc[j] represents a template that would hypothetically appear in fp4

	//First we have to sort the whole structure to make sure we follow the same order as in the sequential code
	//Not fulfilling this may cause to discard some cluster due to another cluster that will also be discarded but hasn't yet

	qsort(bc.bcc, bc.total_bicluster, sizeof(struct bicluster_component), sequential_order);

	// Now we start the discarding process
	for (i = 0; i < bc.total_bicluster; i++){
		for (j = 0; j < i; j++){
			// Skip if clusters were generated by the same process or if the "template" cluster has been discarded
			if (bc.bcc[i].rank == bc.bcc[j].rank || bc.bcc[j].discarded) 
				continue;
			// if (bc.bcc[i].row_pair1 < bc.bcc[j].row_pair1)
			// 	continue;
			// // This condition should never happen as each process will have every combination of row pairs given a row_pair1
			// else if (bc.bcc[i].row_pair1 == bc.bcc[j].row_pair1 && bc.bcc[i].row_pair2 < bc.bcc[j].row_pair2)
			// 	continue;

			// Obtain tsize considering the template (columns) in hexadecimal
			// Since a group of rows that give the same division by 4 would be part of the same
			// digit and the rows in bc.bcc[i].column_number are ordered, we can use this approach
			tsize = 0;
			last_div = -1;
			for (k = 0; k < bc.bcc[i].col_total; k++){
				curr_div = bc.bcc[i].column_number[k] / 4;
				if (curr_div != last_div){
					tsize++;
					last_div = curr_div;
				}
			}

			// Mimic templateFoundInFiles's behaviour
			found_component = 0;
			for (k = 0; k < bc.bcc[i].col_total; k++){
				for (t = 0; t < bc.bcc[j].col_total; t++){
					if (bc.bcc[i].column_number[k] == bc.bcc[j].column_number[t]){
						found_component++;
						break;
					}
				}
			}

			// Lazy disarding
			if (found_component == tsize){
				printf("Discarded template: ");
				for (int z = 0; z < bc.bcc[i].col_total; z++) printf("%d ",bc.bcc[i].column_number[z]);
				printf("\n");
				bc.bcc[i].discarded = 1;
			}

		}
	}
}

// Does what printBicluster did before but now having all the biclusters, instead of going one at a time
void print_final_biclusters(){
	int i,j,id;
	id = 1;
	for (i = 0; i < bc.total_bicluster; i++){
		if (bc.bcc[i].discarded)
			continue;

		fprintf(fp2,"\n bc#: %d",id++);
		fprintf(fp2,"\n c: ");
		for (j = 0; j < bc.bcc[i].col_total; j++){
			fprintf(fp2, "%d ", bc.bcc[i].column_number[j]+1);
		}
		fprintf(fp2, "\n r: ");
		for (j = 0; j < bc.bcc[i].row_total; j++){
			fprintf(fp2, "%d ", bc.bcc[i].row_number[j]+1);
		}
		fprintf(fp2,"\n Total col=%d, Total row=%d",bc.bcc[i].col_total,bc.bcc[i].row_total);
	}
}


//----------------------------
// Converts one int to string format according to a base

char* itxoa(int value, char* buffer, int base)
{
	int n,i,r,j,k=0;
	char buf[100];
   
    // consider the absolute value of the number
    n = abs(value);
    i = 0;
    while (n)
    {
        r = n % base;
		buf[i++] = 48 + r;
        n = n / base;
    }
 
    if (i == 0) {buffer[i++] = '0'; }
	//else
    if (value < 0 && base == 10) { buffer[i++] = '-';   }


	for(j=i-1;j>=0;j--)
		buffer[k++]=buf[j];

	buffer[i] = '\0';

	return buffer;

}

//=================


int main(int argc, char **argv)
{
	// MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0){
		printf("Running pararell RUBIC with %d processes\n", nprocs);
	}

	char s[1000000],s1[10],*s2,fname[200],fname1[400];
	//FILE *fp1,*fp2,*fp3;
	FILE *fp1,*fp3;
	long int i,j,k;   //,**matrix;
	time_t start,end;
	struct timeb st,en;
	int difftm,ccnt;// ccnt and rcnt stands for no. fo clumn and no. of row which are 1
	float a,threshold;
	int max_num_bicluster,extra_col=0;
	char fp4_name[25]; // One file per process for templates

	// Each process gets a separate file for storing its biclusters and templates
	sprintf(fp4_name, "./op_template-%d.txt",rank);
	fp4=fopen(fp4_name,"w+");
	if (fp4 == NULL)
		printf("Problem in opening template file at process %d\n", rank);

	// Files managed by process 0 
	if (rank == 0){
		fp2=fopen(argv[2],"w"); // Last one
	
		if(fp2==NULL)
			puts("Problem in opening output file");
		else
			puts("OK for opening output file");
			
		fprintf(fp2," 1st line shows Bicluster number 'bc#:', 2nd line shows column number 'c:', 3rd line shows row number 'r:' and so on ....");
	
	
		fp3=fopen("./Execution_Result_Pararell.txt","a");
		if(fp3==NULL)
			puts("Problem in opening fp3 output file");
		else
			puts("OK for opening fp3 output file");
	}	

	// Input parameters and file used by all processes
	fp1=fopen(argv[1],"r"); //Last one

	if(fp1==NULL)
		printf("Problem in opening input file at process %d\n", rank);
	else
		printf("OK for opening new_input file by process %d\n",rank);

	mnr=atoi(argv[3]);
	mnc=atoi(argv[4]);

	clst.rsize=0;
	while(fgets(s,1000000,fp1))
		clst.rsize++;

	clst.csize=0;
	for(i=0;i<strlen(s);i++)
		if(s[i]==',')
			clst.csize++;
	clst.csize++;// one commma means 2 column and so so

	extra_col=(clst.csize%4==0)? 0: (4-(clst.csize%4)); // to manage column size as multiple of 4
	// printf("\n clst.rsize=%d, csize=%d",clst.rsize,clst.csize);
	clst.csize+=extra_col; // to manage column size as multiple of 4
	//printf("\n modified csize=%d",clst.csize);
	
	clst.matrix=(int **)calloc(clst.rsize,sizeof(int *));
	for(i=0;i<clst.rsize;i++)
		clst.matrix[i]=(int *)calloc(clst.csize,sizeof(int ));

	i=0;
	rewind(fp1);

	threshold=atof(argv[5]); //Last One
	//printf("\n::: threshold=%f\n",threshold);
	if (rank == 0){
	printf("\n::: mnr=%d\n",mnr);
	printf("\n::: mnc=%d\n",mnc);
	}

	while(fgets(s,1000000,fp1))
	{
			j=0;
			s2=strtok(s,",");
			a=atof(s2);
			
			if(a<threshold)
			{	clst.matrix[i][j]=0; }
			else
			{  clst.matrix[i][j]=1;	}

			while((s2=strtok(NULL,","))!=NULL && ++j<clst.csize)
			{
				a=atof(s2);
				if(a<threshold)
				{	clst.matrix[i][j]=0;}
				else
				{ clst.matrix[i][j]=1;}
			}
			while(j<clst.csize) // to fill the remaining column(extra column) with zero
			{ clst.matrix[i][++j]=0;}
		i++;
	}
	clst.rw=(struct row *)calloc(clst.rsize,sizeof(struct row )); // it will represent each row but with column size=col/4
	clst.row_col=clst.csize/4;
	for(i=0;i<clst.rsize;i++)
	{
		clst.rw[i].rg=(struct row_group *)calloc(clst.row_col,sizeof(struct row_group)); // it will represent each row but with column size=col/4
		// So rg have row_col number of columns Each column represent 4bit decimal number
		for(j=0;j<clst.row_col;j++)
		{
			clst.rw[i].rg[j].value=binaryTodecimal(i,j);//clst.matrix[i][4*j],// indivisual column of rw gets its decimal col value
		}
	}

	// Bicluster bc and sets v initialization
	max_num_bicluster= clst.rsize*(clst.rsize-1)/2 +1;//(clst.rsize*(clst.rsize-1)/2 +1);
	bc.bcc=(struct bicluster_component *)calloc(max_num_bicluster,sizeof(struct bicluster_component ));//old and original
	bc.total_bicluster=0;

	initializeSet();

	// // Sending all necessary information to the remaining processes

	// //mnr and mnc
	// MPI_Bcast(&mnr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// MPI_Bcast(&mnc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// // clst
	// MPI_Bcast(&(clst.csize), 1, MPI_INT, 0, MPI_COMM_WORLD);
	// MPI_Bcast(&(clst.rsize), 1, MPI_INT, 0, MPI_COMM_WORLD);

	// if(rank != 0){
	// 	clst.matrix = (int **)calloc(clst.rsize, sizeof(int *));
	// 	for (int t = 0; t < clst.rsize; t++)
	// 		clst.matrix[t] = (int *)calloc(clst.csize, sizeof(int));
	// }
	// for (int t = 0; t < clst.rsize; t++){
	// 	MPI_Bcast(clst.matrix[t], clst.csize, MPI_INT, 0, MPI_COMM_WORLD);
	// }

	// MPI_Bcast(&(clst.row_col), 1, MPI_INT, 0, MPI_COMM_WORLD);

	// if(rank != 0){
	// 	clst.rw = (struct row *)calloc(clst.rsize, sizeof(struct row));
	// 	for (int t = 0; t < clst.rsize; t++){
	// 		clst.rw[t].rg = (struct row_group *)calloc(clst.row_col, sizeof(struct row_group)); 
	// 	}
	// }

	// // OVERHEAD: for each row group's value (1 int representing 4 bits), we have to broadcast it independentely as it is wrapped in a struct datatype
	// for (int t = 0; t < clst.rsize; t++){
	// 	for (int tt = 0; tt < clst.row_col; tt++){
	// 		MPI_Bcast(&(clst.rw[t].rg[tt].value), 1, MPI_INT, 0, MPI_COMM_WORLD);
	// 	}
	// }

	if (rank == 0){
		puts("\n Calling RUBic");
		ftime(&st);
	}

	// Call to RUBIC
	RUBIC();

	printf("Process %d found %d biclusters\n", rank, num_bicluster);

	// Gather all biclusters
	gather_resuls();

	if (rank == 0){	
		// Discard those that wouldn't appear in the original code
		discard_clusters();
		// Store those clusters in fp2
		print_final_biclusters();
		
		ftime(&en);
		
		// write in output file in pattern 100000010101....

		difftm=(int)(1000*(en.time-st.time)+(en.millitm-st.millitm));
		printf("\n Total time is %u millisecs\nTotal Biclusters %d\n",difftm,num_bicluster);
		fprintf(fp2,"\n\n Total time is %u millisecs",difftm);
		fprintf(fp3,"\n%u,%d",difftm,num_bicluster);

		fclose(fp2);
		fclose(fp3);
	}

	fclose(fp1);
	fclose(fp4);
	remove(fp4_name);

	MPI_Finalize();
}

//-------------------------------
// RUBIC function: executes the biclustering algorithm

void RUBIC()
{
	int i,j,k,kk,l,*t,*tt,nonzero_component_num=0;
	int *temp,temp_sum=0,m=0,n,si=0;
	char s[10],*temp_s,*temp_s1;

	// To control the iterations of the main loop
	int rows_per_proc = (clst.rsize - 1) / 2 / nprocs; // rsize - 1 because the last row doesn't have any pairs
	int start_row = rank * rows_per_proc, end_row_mirror = clst.rsize - 2 - start_row;
	int end_row = (rank+1) * rows_per_proc - 1, start_row_mirror = clst.rsize - 2 - end_row;

	t=(int *)calloc(clst.row_col+1,sizeof(int));
	temp=(int *)calloc(clst.row_col+1,sizeof(int));
	temp_s=(char *)calloc(2*clst.row_col+1,sizeof(char));

	// if (rank != nprocs-1)
	// printf("Process %d: %d-%d and %d-%d\n", rank, start_row, end_row, start_row_mirror, end_row_mirror);

	// if (rank == nprocs-1)
	// printf("Process %d: %d-%d\n", rank,start_row, end_row_mirror);

	i = start_row;
	while( i <= end_row_mirror)
	{
		// Start with the mirror rows (to maintain worload balanced)
		if (rank != nprocs - 1 && i == end_row + 1) // First condition to leave the last process with all the intermediate rows that are left
		{
			i = start_row_mirror;
		}

		m=0; //for each row intialize
		temp_sum=0;
		for(j = i+1; j < clst.rsize; j++)
		{
			printf("Process %d - i,j = %d,%d\n", rank, i+1, j+1);
		    temp_sum=0;
			nonzero_component_num=0;
			for(k=0;k<clst.row_col;k++)// number of column for each row in decimal value not binary
			{
				t[k]=(clst.rw[i].rg[k].value)&(clst.rw[j].rg[k].value); // AND is done to get template in between ith and jth row in decimal value: cluster,row,rowgroup,value
				temp_sum+=t[k];
				if(t[k]!=0)
				{
				  nonzero_component_num++;
				}
			 }

			if(temp_sum!=0)
			{

		//LAST ONE	//if( mncRequirement(t)==1 && !New_templateFound(temp_s))// number of nonzero component is greater than 0; for new bicluster but not all columns are 0
                if( mncRequirement(t)==1 && ! templateFound_inFile(t,nonzero_component_num))//,temp_s))
                {
                    if(search(i,j,t)==0)//;// row, column and template
                     break; // invalid cluster as mnr does not satiesfy
                }
			}//end of if temp_sum
		}// end of for j

	i++;
	}//end of while i
	free(t);
}

//------------------------
// Not used in the code

int New_templateFound(char *t)
{
	char *s,*str;

	s=(char *)calloc(2*clst.row_col+1,sizeof(char));


	rewind(fp4);
	while(fgets(s,strlen(t)+1,fp4)!=NULL)//feof(fp4))//fgets(s,clst.row_col+10,fp4)!=EOF)
	{

		if(strcmp(s,t)==0) //strncmp(s,t,strlen(t)))
		{
		//	puts("MATCHED");
		//	getch();
			free(s);
			return(1);
		}
	}
	free(s);
	return(0);

}

//----------------------------
//----------------------------
// Determine wether the template stored in t was or not previously found searching the correspondent columns in fp4 

int templateFound_inFile(int *t, int tsize)               //,char *temp_s)
{
	int i,j,k,found,found_component=0,sum=0;
	struct bicluster pattern;
	char *s,s1[40],s2[40],*s3,*str;

	pattern.bcc=(struct bicluster_component *)calloc(1,sizeof(struct bicluster_component ));//old and original
	pattern.bcc[0].col_total=0;
	patternConversionDecimalToBinary(t,&pattern);

	//puts("\n Enter in templatefound");

	if(tsize<mnc)//==1 required minimum column is 2 or i.e. <mnc
	 return (1); //invalid template

	//t s=(char *)calloc(2*clst.row_col+1,sizeof(char));
	s=(char *)calloc(4000,sizeof(char));
	//s3=(char *)calloc(4000,sizeof(char));

	for(i=0;i<tsize;i++)
		sum+=t[i];

	if(sum==0)
		return(1);

	rewind(fp4);

	printf("Template of Process %d - t = ",rank);
	for (int i = 0; i < pattern.bcc[0].col_total; i++)
		printf("%d ",pattern.bcc[0].column_number[i]+1);
	printf("\n");

	while(fgets(s,4000,fp4)!=NULL)
	{

		found_component=0;

	//	printf("\n  found_component=%d tsize=%d",found_component,tsize);

		for(k=0;k<pattern.bcc[0].col_total;k++)
		{
			found=0;

			strcpy(s2," ");
			itxoa(pattern.bcc[0].column_number[k]+1,s1,10); //convert column number into base 10 string representation
			strcat(s1," ");
			strcat(s2,s1);


			if(strstr(s,s2)!=NULL)// && strcmp(s,s1)==0)
			{
				found_component++;
			}

		}

		if(found_component==tsize)
		{
			printf("Discarded by template: %s\n",s);
			return(1); //template found
		}

	}

	free(s);
	free(pattern.bcc);
	return(0); // found false  template not found and got new template
}
//----------------------------
// Determines if the min nº of nonzero columns is sattisfied in the template

int mncRequirement(int *t)
{
	int i,j,k,found,found_component;
	struct bicluster pattern;

	pattern.bcc=(struct bicluster_component *)calloc(1,sizeof(struct bicluster_component ));//old and original
	pattern.bcc[0].col_total=0;
	New_patternConversionDecimalToBinary(t,&pattern);

	if(pattern.bcc[0].col_total>=mnc)//==1 required minimum column is 2 or i.e. <mnc
	{
		free(pattern.bcc);
		return (1);
	}
	else
	{
		free(pattern.bcc);
		return(0);
	}
}
//---------------------------
// Not called in the code

int templateFound(int *t)
{
	int i,j,k,found,found_component;
	struct bicluster pattern;

	pattern.bcc=(struct bicluster_component *)calloc(1,sizeof(struct bicluster_component ));//old and original
	pattern.bcc[0].col_total=0;
	patternConversionDecimalToBinary(t,&pattern);

	if(pattern.bcc[0].col_total<mnc)//==1 required minimum column is 2 or i.e. <mnc
	 return (1);


	found=0; // minimum col 2

	bc.total_bicluster=num_bicluster;

	printf("\n bc.total_bicluster=%d",bc.total_bicluster);

	for(i=0;i<bc.total_bicluster;i++)
	{
	//	puts("\n Enter in In templatefound");
		found_component=0;

		for(k=0;k<pattern.bcc[0].col_total;k++)
		{
			found=0;
			for(j=0;j<bc.bcc[i].col_total;j++)
			{

			//if(t[j]==bc.bcc[i].column_number[j])
				if(pattern.bcc[0].column_number[k]==bc.bcc[i].column_number[j])
				{
					found_component++;
					found=1;
		//		printf("\nTemplate Found fn: t[%d]=%d ",j,t[j]);
				}
			
			}//  end for j
		} // end for k


	//t	if(found_component==clst.row_col)// all component of t is found in some previous template stored in templ
	if(found_component==pattern.bcc[0].col_total)// all component of t is found in some previous template stored in templ
		{
			return(1); // found true
		}
	}// end for i

	//}// end if pattern.bcc[0].col_total>1
	free(pattern.bcc);
	return(0); // found false  template not found and got new template
}

//----------------
// Look for rows having the same pattern as r and j

int search(int r,int j,int *t) //
{
	int n,l,m,i,k,*selected_row,selected_row_total=0,*filtered_row,filtered_row_total=0;
	selected_row=(int *)calloc((clst.rsize-2),sizeof(int));

	if(selected_row==NULL)
		puts(" Memory not found for selected_row");

	for(l=1;l<r;l++) // skipping rth row
	{
	//t	if(l!=j)
		selected_row[selected_row_total++]=l; // these row are need to be searched
	}
	for(l=r+1;l<j;l++) // skipping jth row
	{
	//t	if(l!=j)
		selected_row[selected_row_total++]=l; // these row are need to be searched
	}
	for(l=j+1;l<clst.rsize;l++) //As template is between r and j so skip jth row // skipping jth row so need to use	if(l!=j)
	{
		selected_row[selected_row_total++]=l; // these row are need to be searched
	}

	m=0;// 1st column of rg(row_group) of rw(row) of clst(cluster)
	while(m<clst.row_col)
	{
		i=0;
		
		while(i<selected_row_total) //for(l=selected_row[row_num];l<clst.rsize;l++) // going to search new bicluster
		{
			l=selected_row[i];// 1st Column of rg
			if( v[0].set[0]==t[m] || binarySearch(clst.rw[l].rg[m].value,t[m]) )
			{
				selected_row[filtered_row_total++]=l;
			}
			i++;
		}// end while i

		selected_row_total=filtered_row_total;
		filtered_row_total=0;
		m++;// next column of rg(row_group) of rw(row) of clst(cluster)
	}// end while

	if( (selected_row_total+2)<mnr )  // invalid cluster due to less mnr
		return(0);
	else
//t	if((selected_row_total+2)>=mnr) //0
	{
		bc.bcc[bc.total_bicluster].row_number=(int *)calloc(selected_row_total+2,sizeof(int));

		bc.bcc[bc.total_bicluster].row_number[0]=r; // As by for row no. r and j (AND of r and j row gives template) we have worked to find bicluster
		bc.bcc[bc.total_bicluster].row_number[1]=j;

		bc.bcc[bc.total_bicluster].row_pair1 = r;
		bc.bcc[bc.total_bicluster].row_pair2 = j;

		bc.bcc[bc.total_bicluster].rank = rank;
		bc.bcc[bc.total_bicluster].discarded = 0;

		for(i=0;i<selected_row_total;i++)
		{
			bc.bcc[bc.total_bicluster].row_number[i+2]=selected_row[i];
		}

		decimalTobinary_fromTemplate(t);

		bc.bcc[bc.total_bicluster].row_total=selected_row_total+2;

		//if(bc.total_bicluster<1600)
		bc.total_bicluster++; // as one more cluster has been found
		num_bicluster++; //new refferring to total

	}// end of if

	printBicluster();
	
	// free(bc.bcc[--bc.total_bicluster].row_number);
	free(selected_row);
	return (1); //valid cluster

}

//------------------
// Stores bicluster inside file fp5 and fp4
// NOTE: fp5 will store both rows, cols and full information about biclusters. 
// fp4 is used merely to check wether a template o columns has been already included

void printBicluster()
{ 
	int i,j;

        // fprintf(fp5,"\n bc#: %d.%d", rank, num_bicluster); // Include rank in bicluster id
		// fprintf(fp5,"\n c: ");
		fprintf(fp4," ");

		for(j=0;j<bc.bcc[bc.total_bicluster-1].col_total;j++)
		{
			{
				// fprintf(fp5,"%d ",bc.bcc[bc.total_bicluster-1].column_number[j]+1); // print start from 1 not 0
				fprintf(fp4,"%d ",bc.bcc[bc.total_bicluster-1].column_number[j]+1);
			}
		}
		// fprintf(fp5,"\n r: ");
		for(j=0;j<bc.bcc[bc.total_bicluster-1].row_total;j++)
		{
			// fprintf(fp5,"%d ",bc.bcc[bc.total_bicluster-1].row_number[j]+1);
		}// end for j
		// fprintf(fp5,"\n Total col=%d, Total row=%d",bc.bcc[bc.total_bicluster-1].col_total,bc.bcc[bc.total_bicluster-1].row_total);

		fprintf(fp4,"\n");

}

//-----------------
// Checks wether or not value t belongs to the set Vx by searching in the v[x].set array

int binarySearch(int x,int t)
{
	int i,found=0,beg,end,mid;

	if(v[x].total_set<=3) // linear search
	{
		for(i=0;i<v[x].total_set;i++)
		{
			if(v[x].set[i]==t)
			{
				found=1;
				break;
			}
		}// end of for
		return(found);
	}// end of if
	else              // binary search
	{
		beg=0;
		end=v[x].total_set-1;
		mid=(beg+end)/2;
		while(beg<=end)// && v[x].set[mid]!=t)
		{
			if(v[x].set[mid]==t)
				return(1);
			else
			{
				if(t<v[x].set[mid])
				 end=mid-1;
				else
				 beg=mid+1;
				 mid=(beg+end)/2;
			}
		}

	return(0);
	}// end of else

}


//--------------------------
// Similar to PatternConversionDecimalToBinary

void New_patternConversionDecimalToBinary(int *t,struct bicluster *pattern) // pattern.col store only those column number which corresponds to 1.
{																		// pattern.col_total will store total number of column having 1 not 0.
	int i,j,remainder,tt,cnt,*temp,nonzero=0;

	temp=(int *)calloc(clst.csize,sizeof(int));

	for(i=0;i<clst.row_col;i++)
	{
		tt=t[i];
		for(cnt=0;cnt<4;cnt++)
		{
			remainder=tt%2;
			tt/=2;
			temp[(i+1)*4-(1+cnt)]=remainder;
			if(remainder!=0)
				nonzero++;
		}
	}

	pattern->bcc[0].col_total=nonzero;
	free(temp);
}

//-----------------
//Takes the template in decimal, convertis it to binary and stores it in the bicluster set

void patternConversionDecimalToBinary(int *t,struct bicluster *pattern) // pattern.col store only those column number which corresponds to 1.
{																		// pattern.col_total will store total number of column having 1 not 0.
	int i,j,remainder,tt,cnt,*temp,nonzero=0;

	temp=(int *)calloc(clst.csize,sizeof(int));

	for(i=0;i<clst.row_col;i++)
	{
		tt=t[i];
		for(cnt=0;cnt<4;cnt++)
		{
			remainder=tt%2;
			tt/=2;
			temp[(i+1)*4-(1+cnt)]=remainder;
			if(remainder!=0)
				nonzero++;
		}
	}

	pattern->bcc[0].column_number=(int *)calloc(nonzero,sizeof(int));

	pattern->bcc[0].col_total=nonzero;
	i=0;
	for(j=0;j<clst.csize;j++) // will say which column are 1
		{
			if(temp[j]==1)
			{
				pattern->bcc[0].column_number[i++]=j; // storing col number
			}
		}
	free(temp);
}

//--------------------
// Converts the values of a template in binary format and stores the original column numbers of the binary matrix in the bicluster set

void decimalTobinary_fromTemplate(int *t)
{
	int i,j,remainder,tt,cnt,*temp,nonzero=0;
	temp=(int *)calloc(clst.csize,sizeof(int));

	for(i=0;i<clst.row_col;i++)
	{
		tt=t[i];
		for(cnt=0;cnt<4;cnt++)
		{
			remainder=tt%2;
			tt/=2;
			temp[(i+1)*4-(1+cnt)]=remainder;
			if(remainder!=0)
				nonzero++;
		}
	}
	bc.bcc[bc.total_bicluster].column_number=(int *)calloc(nonzero,sizeof(int));
	bc.bcc[bc.total_bicluster].col_total=nonzero;
	i=0;
	for(j=0;j<clst.csize;j++) // will say which column are 1
		{
			if(temp[j]==1)
			{
				bc.bcc[bc.total_bicluster].column_number[i++]=j; // storing col number
			}
		}
	free(temp);
}

//---------------
//Not used in the code

void decimalTobinary(struct bicluster_component *b)
{
	int i,j,remainder,tt,cnt,*t;

	t=(int *)calloc(clst.row_col,sizeof(int));

	for(i=0;i<clst.row_col;i++)
		t[i]=b->column_number[i];

	for(i=0;i<clst.row_col;i++)
	{
		tt=t[i];
		for(cnt=0;cnt<4;cnt++)
		{
			remainder=tt%2;
			tt/=2;
			b->column_number[(i+1)*4-(1+cnt)]=remainder;
		}
	}
	free(t);

}

//------------------
// Used in main to represent a sequence of 4 consecutive bits of matrix as a decimal value inside of the row_groups

int binaryTodecimal(int r, int j)
{
	int bin[4],dec=0,i;

	for(i=0;i<4;i++)
	{
		dec+=pow(2,3-i)*clst.matrix[r][4*j+i];
	}

	return(dec);
}

//---------------------
// Initialize the 16 sets defined in the documentation as Di

void initializeSet()
{
	v[0].total_set=1;
	v[0].set=(int *)calloc( v[0].total_set,sizeof(int));
	v[0].set[0]=0;

	v[1].total_set=1;
	v[1].set=(int *)calloc( v[1].total_set,sizeof(int));
	v[1].set[0]=1;


	v[2].total_set=1;
	v[2].set=(int *)calloc( v[2].total_set,sizeof(int));
	v[2].set[0]=2;

	v[3].total_set=3;
	v[3].set=(int *)calloc( v[3].total_set,sizeof(int));
	v[3].set[0]=1;
	v[3].set[1]=2;
	v[3].set[2]=3;


	v[4].total_set=1;
	v[4].set=(int *)calloc( v[4].total_set,sizeof(int));
	v[4].set[0]=4;


	v[5].total_set=3;
	v[5].set=(int *)calloc( v[5].total_set,sizeof(int));
	v[5].set[0]=1;
	v[5].set[1]=4;
	v[5].set[2]=5;


	v[6].total_set=3;
	v[6].set=(int *)calloc( v[6].total_set,sizeof(int));
	v[6].set[0]=2;
	v[6].set[1]=4;
	v[6].set[2]=6;


	v[7].total_set=7;
	v[7].set=(int *)calloc( v[7].total_set,sizeof(int));
	v[7].set[0]=1;
	v[7].set[1]=2;
	v[7].set[2]=3;
	v[7].set[3]=4;
	v[7].set[4]=5;
	v[7].set[5]=6;
	v[7].set[6]=7;


	v[8].total_set=1;
	v[8].set=(int *)calloc( v[8].total_set,sizeof(int));
	v[8].set[0]=8;


	v[9].total_set=3;
	v[9].set=(int *)calloc( v[9].total_set,sizeof(int));
	v[9].set[0]=1;
	v[9].set[1]=8;
	v[9].set[2]=9;


	v[10].total_set=3;
	v[10].set=(int *)calloc( v[10].total_set,sizeof(int));
	v[10].set[0]=2;
	v[10].set[1]=8;
	v[10].set[2]=10;


	v[11].total_set=7;
	v[11].set=(int *)calloc( v[11].total_set,sizeof(int));
	v[11].set[0]=1;
	v[11].set[1]=2;
	v[11].set[2]=3;
	v[11].set[3]=8;
	v[11].set[4]=9;
	v[11].set[5]=10;
	v[11].set[6]=11;


	v[12].total_set=3;
	v[12].set=(int *)calloc( v[12].total_set,sizeof(int));
	v[12].set[0]=4;
	v[12].set[1]=8;
	v[12].set[2]=12;



	v[13].total_set=7;
	v[13].set=(int *)calloc( v[13].total_set,sizeof(int));
	v[13].set[0]=1;
	v[13].set[1]=4;
	v[13].set[2]=5;
	v[13].set[3]=8;
	v[13].set[4]=9;
	v[13].set[5]=12;
	v[13].set[6]=13;


	v[14].total_set=7;
	v[14].set=(int *)calloc( v[14].total_set,sizeof(int));
	v[14].set[0]=2;
	v[14].set[1]=4;
	v[14].set[2]=6;
	v[14].set[3]=8;
	v[14].set[4]=10;
	v[14].set[5]=12;
	v[14].set[6]=14;



	v[15].total_set=15;
	v[15].set=(int *)calloc(v[15].total_set,sizeof(int));
	v[15].set[0]=1;
	v[15].set[1]=2;
	v[15].set[2]=3;
	v[15].set[3]=4;
	v[15].set[4]=5;
	v[15].set[5]=6;
	v[15].set[6]=7;
	v[15].set[7]=8;
	v[15].set[8]=9;
	v[15].set[9]=10;
	v[15].set[10]=11;
	v[15].set[11]=12;
	v[15].set[12]=13;
	v[15].set[13]=14;
	v[15].set[14]=15;


}
