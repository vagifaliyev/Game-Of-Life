// -------- PACKAGES --------
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <math.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstdint>

#define periodic true // <---------------------------------------- BOUNDARY CONDITION can be changed here

using namespace std; 

// -------- GLOBAL VARIABLES --------

int id; // processor identification number
int p; // total no of processors.
int sp_rows; // no of total rows in a specific processor
int sp_cols; // no of total columns in a specific processor

//decleration of domains
bool **domain;
bool **new_domain;

// g stands for Global
int g_imax = 5000;	// <---------------------------------------- GLOBAL I_MAX can be changed here
int g_jmax = 5000;	// <---------------------------------------- GLOBAL J_MAX can be changed here

int l_imax, l_jmax; // l stand for Local

int no_steps = 10; // <-------------------------------------------- NUMBER OF STEPS can be changed here 

int d_row; // the row of a processor in the domain
int d_col; // the column of a processor in the domain


// -------- FIND DIMENSIONS --------

// get dimensions for processors 
// creates domains for processors
// default way to do it but could be done better
void find_dimensions()
{
	int min_gap = g_imax;

	if (g_jmax > g_imax)
		min_gap = g_jmax;

	//loop until the minimum total processors can be divided
	for (int i = 1; i <= p; i++)
	{
		if (p%i == 0)
		{
			int gap = abs(g_jmax / (p / i) - g_imax / i);

			if (gap < min_gap)
			{
				min_gap = gap;

				// equals iteration number
				sp_rows = i; 
				//number of colums equal total number of processor divided by iteration number
				sp_cols = p / i;
			}
		}
	}

	//print out the dimensions of the domain created 
	if (id == 0)
	{
		cout << "----------------------------" << endl;
		cout << "Divide " << p << " into " << sp_rows << " by " << sp_cols << " grid" << endl;
		cout << "----------------------------" << endl;
	}
}

// -------- MESH GRID -------
//function to determine the cell size(x,y) in the local domain of a processor
void cells_for_grids()
{
	d_row = id / sp_cols;
	d_col = id % sp_cols;

	int rows_remain = g_imax;
	for (int i = 0; i <= d_row; i++){
		l_imax = rows_remain / (sp_rows - i);
		rows_remain -= l_imax;
	}

	int cols_remain = g_jmax;
	for (int j = 0; j <= d_col; j++){
		l_jmax = cols_remain / (sp_cols - j);
		cols_remain -= l_jmax;
	}
}

//function to convert processors id to it's index on the domain
//example id 0's index will be (0,0)
void id_to_index(int id, int &id_row, int &id_column)
{
	//get the column index for the processor id as specified in the argument
	id_column = id % sp_cols;
	//get the row idex for the processor id as specified in the argument
	id_row = id / sp_cols;
}

//function to determine processor's id from its index.
//example index (0,0) will be id 0's 
int id_from_index(int id_row, int id_column)
{
	if (periodic) // PERIODIC domain
	{
		//return processor's id
		return ((id_row + sp_rows) % sp_rows) * sp_cols + (id_column + sp_cols) % sp_cols;
	}
	else if (!periodic)  // NON-PERIODIC domain
	{
		//check if the processor's ROW id is out of bounds
		if (id_row >= sp_rows || id_row < 0)
			return -1;
		//check if the processor's COLUMN id is out of bounds
		if (id_column >= sp_cols || id_column < 0)
			return -1;
		//return processor's id
		return id_row * sp_cols + id_column;
	}
}

void make_comm()
{
	// --------------- DATA TYPES ---------------
	// class implementation fails, using this method instead 

	vector<vector<int>> block_length(16);
	vector<vector<MPI_Aint>> addresses(16);
	vector<vector<MPI_Datatype>> typelist(16);
	vector<MPI_Datatype> temp_data(16);

	// --------- TOP ---------

	//send
	block_length[1].push_back(l_jmax);
	MPI_Aint temp1;
	MPI_Get_address(&domain[1][1], &temp1);
	addresses[1].push_back(temp1);
	typelist[1].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[1].size(), &block_length[1][0], &addresses[1][0], &typelist[1][0], &temp_data[1]);
	MPI_Type_commit(&temp_data[1]);


	//recv
	block_length[9].push_back(l_jmax);
	MPI_Aint temp9;
	MPI_Get_address(&domain[0][1], &temp9);
	addresses[9].push_back(temp9);
	typelist[9].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[9].size(), &block_length[9][0], &addresses[9][0], &typelist[9][0], &temp_data[9]);
	MPI_Type_commit(&temp_data[9]);


	// --------- BOTTOM ---------
	
	//send
	block_length[6].push_back(l_jmax);
	MPI_Aint temp6;
	MPI_Get_address(&domain[l_imax][1], &temp6);
	addresses[6].push_back(temp6);
	typelist[6].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[6].size(), &block_length[6][0], &addresses[6][0], &typelist[6][0], &temp_data[6]);
	MPI_Type_commit(&temp_data[6]);

	//recv
	block_length[14].push_back(l_jmax);
	MPI_Aint temp14;
	MPI_Get_address(&domain[l_imax+1][1], &temp14);
	addresses[14].push_back(temp14);
	typelist[14].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[14].size(), &block_length[14][0], &addresses[14][0], &typelist[14][0], &temp_data[14]);
	MPI_Type_commit(&temp_data[14]);

	// --------- RIGHT ---------

	//SEND
	for (int i = 1; i <=l_imax; i++)
	{
		block_length[4].push_back(1);
		MPI_Aint temp4;
		MPI_Get_address(&domain[i][l_jmax], &temp4);
		addresses[4].push_back(temp4);
		typelist[4].push_back(MPI_C_BOOL);
	}
	
	MPI_Type_create_struct(block_length[4].size(), &block_length[4][0], &addresses[4][0], &typelist[4][0], &temp_data[4]);
	MPI_Type_commit(&temp_data[4]);


	//RECV
	for (int i = 1; i <=l_imax; i++)
	{
		block_length[12].push_back(1);
		MPI_Aint temp12;
		MPI_Get_address(&domain[i][l_jmax+1], &temp12);
		addresses[12].push_back(temp12);
		typelist[12].push_back(MPI_C_BOOL);
	}
	
	MPI_Type_create_struct(block_length[12].size(), &block_length[12][0], &addresses[12][0], &typelist[12][0], &temp_data[12]);
	MPI_Type_commit(&temp_data[12]);

	// --------- LEFT ---------

	//SEND
	for (int i = 1; i <=l_imax; i++)
	{
		block_length[3].push_back(1);
		MPI_Aint temp3;
		MPI_Get_address(&domain[i][1], &temp3);
		addresses[3].push_back(temp3);
		typelist[3].push_back(MPI_C_BOOL);
	}
	
	MPI_Type_create_struct(block_length[3].size(), &block_length[3][0], &addresses[3][0], &typelist[3][0], &temp_data[3]);
	MPI_Type_commit(&temp_data[3]);

	//RECV
	for (int i = 1; i <=l_imax; i++)
	{
		
		block_length[11].push_back(1);
		MPI_Aint temp11;
		MPI_Get_address(&domain[i][0], &temp11);
		addresses[11].push_back(temp11);
		typelist[11].push_back(MPI_C_BOOL);
	}
	
	MPI_Type_create_struct(block_length[11].size(), &block_length[11][0], &addresses[11][0], &typelist[11][0], &temp_data[11]);
	MPI_Type_commit(&temp_data[11]);


	// --------- TOP_RIGHT ---------

	//SEND
	block_length[2].push_back(1);
	MPI_Aint temp2;
	MPI_Get_address(&domain[1][l_jmax], &temp2);
	addresses[2].push_back(temp2);
	typelist[2].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[2].size(), &block_length[2][0], &addresses[2][0], &typelist[2][0], &temp_data[2]);
	MPI_Type_commit(&temp_data[2]);

	//RECV	
	block_length[10].push_back(1);
	MPI_Aint temp10;
	MPI_Get_address(&domain[0][l_jmax+1], &temp10);
	addresses[10].push_back(temp10);
	typelist[10].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[10].size(), &block_length[10][0], &addresses[10][0], &typelist[10][0], &temp_data[10]);
	MPI_Type_commit(&temp_data[10]);

	// --------- TOP_LEFT ---------

	//SEND
	block_length[0].push_back(1);
	MPI_Aint temp0;
	MPI_Get_address(&domain[1][1], &temp0);
	addresses[0].push_back(temp0);
	typelist[0].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[0].size(), &block_length[0][0], &addresses[0][0], &typelist[0][0], &temp_data[0]);
	MPI_Type_commit(&temp_data[0]);

	//RECV	
	block_length[8].push_back(1);
	MPI_Aint temp8;
	MPI_Get_address(&domain[0][0], &temp8);
	addresses[8].push_back(temp8);
	typelist[8].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[8].size(), &block_length[8][0], &addresses[8][0], &typelist[8][0], &temp_data[8]);
	MPI_Type_commit(&temp_data[8]);

	// --------- BOTTOM_LEFT ---------

	//SEND
	block_length[5].push_back(1);
	MPI_Aint temp5;
	MPI_Get_address(&domain[l_imax][1], &temp5);
	addresses[5].push_back(temp5);
	typelist[5].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[5].size(), &block_length[5][0], &addresses[5][0], &typelist[5][0], &temp_data[5]);
	MPI_Type_commit(&temp_data[5]);


	//RECV	
	block_length[13].push_back(1);
	MPI_Aint temp13;
	MPI_Get_address(&domain[l_imax+1][0], &temp13);
	addresses[13].push_back(temp13);
	typelist[13].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[13].size(), &block_length[13][0], &addresses[13][0], &typelist[13][0], &temp_data[13]);
	MPI_Type_commit(&temp_data[13]);


	// --------- BOTTOM_RIGHT ---------

	//SEND
	block_length[7].push_back(1);
	MPI_Aint temp7;
	MPI_Get_address(&domain[l_imax][l_jmax], &temp7);
	addresses[7].push_back(temp7);
	typelist[7].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[7].size(), &block_length[7][0], &addresses[7][0], &typelist[7][0], &temp_data[7]);
	MPI_Type_commit(&temp_data[7]);


	//RECV	
	block_length[15].push_back(1);
	MPI_Aint temp15;
	MPI_Get_address(&domain[l_imax+1][l_jmax+1], &temp15);
	addresses[15].push_back(temp15);
	typelist[15].push_back(MPI_C_BOOL);
	
	MPI_Type_create_struct(block_length[15].size(), &block_length[15][0], &addresses[15][0], &typelist[15][0], &temp_data[15]);
	MPI_Type_commit(&temp_data[15]);
	
	// --------------- DATA TYPES FINISHED ---------------

	int cnt_type_send = 0;

	MPI_Request* request = new MPI_Request[8 * 2];

	int cnt = 0;
	int tag_num = 0;
	int iter = 0; //for skipping self 

	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			int com_i = d_row + i;
			int com_j = d_col + j;
			int com_id = id_from_index(com_i, com_j);
			//cout << "com_i-> " << com_i << "  com_j-> " << com_j << "  com_id-> " << com_id << endl;
			if (iter != 4){
				if (com_id >= 0 && com_id < p){
					MPI_Isend(MPI_BOTTOM, 1, temp_data[cnt_type_send], com_id, tag_num, MPI_COMM_WORLD, &request[cnt * 2]);
					MPI_Irecv(MPI_BOTTOM, 1, temp_data[cnt_type_send + 8], com_id, (8 - tag_num), MPI_COMM_WORLD, &request[cnt * 2 + 1]);
					cnt++;
					// if (id==0){cout << id << "====>" << com_id << " |" << tag_num;
					// cout << "--->"<<8-tag_num << endl;
					// }
					// if (id==2){cout << "pr: "<< id << "====>" << com_id << "|" << tag_num;
					// cout << "++++>"<<8-tag_num << endl;
					// }				
					}
				cnt_type_send++;
			}
			tag_num++;
			iter++;
		}
	}
	MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);
	// MPI_Barrier(MPI_COMM_WORLD);


	for (int i = 0; i < 16; i++){
		MPI_Type_free(&temp_data[i]);}
}


// -------- FIND NO of ALIVE NEIGHBORS --------
//function to count number of neighbours alive
int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++) //loop around all the 8 (eight) neighbours
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0) //ignore itself for counting
			{
				ix = i + ii; //setting the row to be checked
				jx = j + jj; //setting the column to be checked
				//if alive
				if (domain[ix][jx]) cnt++; //increase count
			}
	return cnt; //return total number of neighbours alive
}

// -------- GAME OF LIFE RULES --------
//iteration function - checker 
void do_iteration(void)
{
	//loop through the local domain of the processor
	for (int i = 1; i <= l_imax; i++)
		for (int j = 1; j <= l_jmax; j++)
		{
			//let the new status (dead/alive) of the cell equal previous status
			new_domain[i][j] = domain[i][j];
			//get the number of neighbours of the alive cell
			int num_n = num_neighbours(i, j);
			if (domain[i][j])
			{
				if (num_n != 2 && num_n != 3) // lives if 2/3 neighbours
					// if the cell has fewer than 2 neighbours - 4 or more neighbours, set dead
					new_domain[i][j] = false; 

			}
			else if (num_n == 3) new_domain[i][j] = true; // comes to live if 3 neighbours
		}
	//update the domain before the next iteration.
	bool **temp;
	temp = domain;
	domain = new_domain;
	new_domain = temp;
}

// -------- WRITE OUT --------
//function to output the results from each processor after every iteration
void grid_to_file(int it)
{
	stringstream fname;
	fstream f1;

	fname << "it_" << it << "_row_"<< d_row << "_col_" << d_col << ".txt";

	//open file
	f1.open(fname.str().c_str(), ios_base::out);

	//loop through the processor's local domain and print
	for (int i = 1; i <= l_imax; i++)
	{
		for (int j = 1; j <= l_jmax; j++)
			f1 << domain[i][j] << "  ";
		f1 << endl;
	}

	f1.close(); //close file
}


int main(int argc, char *argv[])
{
	auto start = chrono::high_resolution_clock::now();

	//setting up MPI communications
	MPI_Init(&argc, &argv);
	//Reads the number of the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	//Reads the total number of processes that have been assigned
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//Generation of random values
	srand(time(NULL) + id * 10000);

	//find the dimensions the total number of processors will be divided into
	find_dimensions();

	//calculate the cell size(x, y) in the local domain of each processor
	cells_for_grids();

	//create the indices of the processors in the global domain from their id
	id_to_index(id, d_row, d_col);

	//allocate memory for the local domain i_max of each processor - including buffers
	domain = new bool*[l_imax + 2];
	//allocate memory for the local new_domain i_max of each processor - including buffers
	//this is used for iteration.
	new_domain = new bool*[l_imax + 2];

	//allocate memory for the local domain j_max of each processor - including buffers
	for (int i = 0; i < l_imax + 2; i++)
		domain[i] = new bool[l_jmax + 2];

	//allocate memory for the local new_domain j_max of each processor - including buffers
	//this is needed for iteration.
	for (int i = 0; i < l_imax + 2; i++)
		new_domain[i] = new bool[l_jmax + 2];

	// set the domain with random values of 0 and 1 only. 
	// 1 means alive (true) <-> 0 means dead (false)
	for (int i = 1; i < l_imax + 1; i++){
		for (int j = 1; j < l_jmax + 1; j++)
		{
			//looping around the loacal domain of each processor - excluding the buffer
			domain[i][j] = (rand() % 2);
			new_domain[i][j] = 0;
		}
	}

	for (int n = 0; n < no_steps; n++)
	{
		//setting up communication with neighbouring processors 
		//peer to peer communication

		make_comm();

		//print out the number of iterations and corresponding cores
		cout << " iteration: " << n << " core: "<< id << endl;
		do_iteration(); //doing iteration to implement the game of live 

		grid_to_file(n); //output file for each iteration/core to a txt file
	}
	
	//delete pointers to pointers for domain
	for (int i = 0; i < l_imax + 2; i++)
		delete[] domain[i];
	delete[] domain;

	//delete pointers to pointers for new_domain
	for (int i = 0; i < l_imax + 2; i++)
		delete[] new_domain[i];
	delete[] new_domain;

	auto finish = chrono::high_resolution_clock::now();

	if (id == 0 ){
		cout << "No of Processors --->  " << p << endl;
		cout << "Dimensions: RowxColoumn --->  " << g_imax << " x " << g_jmax << endl;
		std::chrono::duration<double> time_taken = finish - start;
	 	cout << setprecision(5);
		cout << "Time Taken --->" << time_taken.count() << " s" << endl;
	}

	MPI_Finalize();
}