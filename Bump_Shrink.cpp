// Copyright 2017 A.C. Balazs Group @ University of Pittsburgh.  All Rights Reserved.
// Shrink gel
// First version developed by Olga Kuksenok, olga.kuksenok@gmail.com
// Previous version edited by Awaneesh Singh, awaneesh11@gmail.com
// Thank Santidan Biswas for meaningful discussion on looping over elements not nodes
// Previous version edited by Tao Zhang, zhangtao.scholar@gmail.com
// Modified by Ya Liu for gel attach to rigid bump
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>

using namespace std;

#define AddMobileCiliaOn 1
#define Bcilia 3.75
#define Kbend  37.5

// define SP&LIGHT properties
#define SPUNIFORM 0
#define Kdark  0.05
#define ILIGHT 0.5     // (Kdark*10),
#define SPTOT  0.2     // total
#define kappa  0.08250 // Chemomechanical interaction between SP in spiro form and gel, alpha in paper
#define TimeLightON 20000.0

// gel paramerers
#define c0 0.00072  // crosslink density in the undeformed state
//#define phi0 0.3  // initial polymer volume fraction 0.114
#define chi1 0.518  // gel-solvent interaction parameter

//#define deltaX (1.0)

#define L0tmp 300.
#define UseCosHarmonicBend 1

#define MobilityonTemperatureON 1

//bump parameters: 40*10*16
double Lx = 80, Ly = 10, Lz = 32;
double bump_top_thick = 2.;
double bump_height = Lz - bump_top_thick;
double bump_radius = ((0.5*Lx)*(0.5*Lx) + bump_height*bump_height)/(2*bump_height);
double bump_center_x = 0.5*Lx;
double bump_center_y = bump_height-bump_radius;
double deltaX = 1.0;
double PHIIni = 0.025;//0.0503; //20-32
//double PHIIni = 0.124282; //32-20
double phi0 = 0.025;//0.0503;
double LamdaIni = pow((phi0 / PHIIni), (1.0 / 3.0));
double unit_space = deltaX * LamdaIni;

// DEFINITIONS: End ////////////////////////////////////////

class Node;
class Element;

class Node {
public:
	static Element ** elements_array;
	double position[3];
	double FSpring[3];
	double FPressure[3];
	double FCilia[3];
	double FCiliaBend[3];
	double velocity[3];
	double mobility;
	long int elements_index[8];
	int num_elements;
	long int on_fibers;
	int Pos_Node; // = 0:bump; = 1:edge gel; =2:center gel
	Node();
	double compute_phi(); //AverageBulk
};

class Element {
public:
	static Node ** nodes_array;
	long int nodes_index[8];
	double phi;
	double volume;
	double vSP;
	double vSPTOT;
	double P0;
	Element();
	void update_volume();
	void update_nodes_FSpring();
	void update_nodes_FPressure();
};

Node::Node() {
	for (int i = 0; i < 3; i++) {
		position[i] = 0.;
		FSpring[3] = 0.;
		FPressure[3] = 0.;
		FCilia[i] = 0.;
		FCiliaBend[i] = 0.;
		velocity[i] = 0.;
	}
	mobility = 0.;
	num_elements = 0;
	on_fibers = -1;
	for (int i = 0; i < 8; i++)
		elements_index[i] = -1;
}

double Node::compute_phi() {
	double tmp = 0.0;
	int i;
	long int element_index;

	for (i = 0; i < num_elements; i++) {
		element_index = elements_index[i];
		tmp += elements_array[element_index]->phi;
	}

	tmp /= num_elements;

	return tmp;
}

Element::Element() {
	phi = 0.;
	P0 = 0.;
	vSP = 0.;
	vSPTOT = 0.0;
}

int    Initialize_All();
void   isBump();//only apply at beginning
//void   Read_XYZ();
void   Update_XYZ();
void   CalculateVolumes();
void   CalculateMobilities();
void   CalculatePressure(double);
void   CalculateVelocities();
void dump_configuration_vtk(int);
void log_time_phi();
void log_node_coordinates(long int node_index);
void log_node_velocities(long int node_index);
void log_node_forces(long int node_index);
int  load_conf(char *filename);

//double PHIIni;   // initial phi
double   t_start, t_end;
double   h;
double   temperature;
int   Nfiber;
int **fibers;
long int **nodes_on_fibers;
int *fibers_length;
double chi0;
double Vol0;
double Vol0_phi0;
double c06;
double fiber_equi_length;
double L0;
double    LightSpace;
long int simulation_step;
long int count_dump, count_dump_0;
long int count_log, count_log_0;
int dump_file_type;
Node ** nodes;
Element ** elements;
long int N_nodes;
long int N_elements;

Node ** Element::nodes_array = NULL;
Element ** Node::elements_array = NULL;

int main(int argc, char *argv[]) {
	long int i, j, k;
	double   simulation_time;

	Initialize_All();
	cout<<"check if it is bump"<<endl;
	//isBump();
	// Done ///////////////////////////////////////////////////////////////

	// Start simulation for each time step ////////////////////////////////
	simulation_step = 0;

	count_dump = count_dump_0;
	count_log_0 = 1000;
	count_log = count_log_0;

	simulation_time = t_start;
	cout << "Start Simulation" << endl;
	int t_label = 0;
	while (simulation_time <= t_end) {
		if (count_dump == count_dump_0) {
			count_dump = 0;
			dump_configuration_vtk(t_label);
			t_label++;
			cout << simulation_time << " done..." << endl;
		}

		count_dump++;

		//cout << "t: " << simulation_time << endl;
		CalculateVolumes();  // Calculates volumes of all elements and update phi
		//cout << "checkpoint 1" << endl;
		CalculatePressure(simulation_time);
		//cout << "checkpoint 2" << endl;
		CalculateMobilities();
		//cout << "checkpoint 3" << endl;
		CalculateVelocities();  // Calculates nodal velocities
		//cout << "checkpoint 4" << endl;

		if (count_log == count_log_0) {
			count_log = 0;
			// log_time_phi();
			log_node_coordinates(22);
			log_node_velocities(22);
			log_node_forces(22);
			log_node_coordinates(652);
			log_node_velocities(652);
			log_node_forces(652);
			cout << simulation_time << " log data" << endl;
		}

		count_log++;

		Update_XYZ();  // Updates nodal positions
		//cout << "checkpoint 5" << endl;

		simulation_step++;
		simulation_time = simulation_step * h;
	}
	cout << "Simulation Done" << endl;
	// Done ///////////////////////////////////////////////////////////////////

	// free memory space /////////////////////////////////////////////////////
	if (Nfiber > 0) {
		for (i = 0; i < Nfiber; i++) {
			free(nodes_on_fibers[i]);
		}

		free(nodes_on_fibers);
		free(fibers_length);
	}
	cout << "Free Space Done" << endl;
	// Done //////////////////////////////////////////////////////////////////

	return 0;
}

int Initialize_All() {
	long int i, j, k;
	int l;
	long int node_index, element_index;

	Nfiber = 0;

	load_conf((char *)"swell.conf");
	cout << "loaded swelling conf file" << endl;

	if (fabs(temperature - 20.0) < 1e-6){
		PHIIni = 0.124282; // 32 -> 20
	}
	else if (fabs(temperature - 32.0) < 1e-6) {
		PHIIni = 0.025;//0.0503; // 20 -> 32
	}
	else {
		cout << "Wrong temperature." << endl;
		exit(1);
	}


	c06 = c0 * deltaX / 12.0;

	Vol0 = (deltaX * deltaX * deltaX);

	Vol0_phi0 = Vol0 * phi0;

	chi0 = (30.614891304347623 + 3.4181159420289857 * temperature) / (273.15 + temperature);

	fiber_equi_length = deltaX * (pow((phi0 / PHIIni), (1. / 3.)));

	L0 = (sqrt(phi0) * L0tmp) / (deltaX * deltaX * deltaX);

	ifstream vtkfile("fullsample.gel.vtk");


	if (!vtkfile.is_open())
	{
		cout << "Error opening sample vtk file" << endl;
		exit(1);
	}

	char buffer[256];
	long int tmp_long;

	vtkfile.getline(buffer, 100);
	vtkfile.getline(buffer, 100);
	vtkfile.getline(buffer, 100);
	vtkfile.getline(buffer, 100);
	vtkfile >> buffer;
	vtkfile >> N_nodes;
	cout << "N_nodes: " << N_nodes << endl;
	vtkfile >> buffer;

	// initialize nodes objects
	nodes = (Node **)malloc(sizeof(Node *)*N_nodes);
	for (node_index = 0; node_index < N_nodes; node_index++) {
		nodes[node_index] = new Node;
	}

	Element::nodes_array = nodes;

	// read positions of nodes
	/*
	for (i = 0; i < N_nodes; i++)
		for (j = 0; j < 3; j++) {
			vtkfile >> nodes[i]->position[j];
	}
	*/
	for (i = 0; i < N_nodes; i++){
		vtkfile >> nodes[i]->position[0]>> nodes[i]->position[1]
		>>nodes[i]->position[2]>>nodes[i]->Pos_Node;
	}
	vtkfile.getline(buffer, 100);
	vtkfile >> buffer;
	vtkfile >> N_elements;
	cout << "N_elements: " << N_elements << endl;
	vtkfile >> tmp_long;

	// initialize elements objects
	elements = (Element **)malloc(sizeof(Element *)*N_elements);
	for (element_index = 0; element_index < N_elements; element_index++) {
		elements[element_index] = new Element();
	}

	Node::elements_array = elements;

	/////////////// set up relations between nodes and elements /////////////
	for (element_index = 0; element_index < N_elements; element_index++) {
		vtkfile >> tmp_long;
		//cout << tmp_long << " ";
		for (i = 0; i < tmp_long; i++) {
			vtkfile >> elements[element_index]->nodes_index[i];
			//cout << elements[element_index]->nodes_index[i] << " ";
		}
		//cout << endl;
	}

	for (element_index = 0; element_index < N_elements; element_index++) {
		for (i = 0; i < 8; i++) {
			node_index = elements[element_index]->nodes_index[i];
			nodes[node_index]->elements_index[nodes[node_index]->num_elements] = element_index;
			nodes[node_index]->num_elements = nodes[node_index]->num_elements + 1;
		}
	}
	////////////////////////////////////////////////////////////////////////////////

	vtkfile.close();

	for (element_index = 0; element_index < N_elements; element_index++) {
		elements[element_index]->update_volume();
		elements[element_index]->phi = Vol0_phi0 / elements[element_index]->volume;
		if (SPUNIFORM == 1)
			elements[element_index]->vSPTOT = SPTOT;
	}

	// Initialize fibers ////////////////////////////////////////////////////
	ifstream fibers_vtkfile("sample.fibers.vtk");

	if (!fibers_vtkfile.is_open())
	{
		cout << "Error opening sample fibers vtk file" << endl;
		//exit(1);
	}

	fibers_vtkfile.getline(buffer, 100);
	fibers_vtkfile.getline(buffer, 100);
	fibers_vtkfile.getline(buffer, 100);
	fibers_vtkfile.getline(buffer, 100);
	fibers_vtkfile >> buffer;
	fibers_vtkfile >> tmp_long;
	fibers_vtkfile >> buffer;

	for (i = 0; i < tmp_long; i++) {
		fibers_vtkfile.getline(buffer, 100);
	}
	fibers_vtkfile.getline(buffer, 100);

	fibers_vtkfile >> buffer;
	fibers_vtkfile >> Nfiber;
	fibers_vtkfile >> tmp_long;
//  Ya's change: no fiber
    Nfiber = 0;
	if (Nfiber > 0) {
		nodes_on_fibers = (long int **)malloc(sizeof(long int *) * Nfiber);
		fibers_length = (int *)malloc(sizeof(int) * Nfiber);

		for (int l = 0; l < Nfiber; l++) {
			fibers_vtkfile >> fibers_length[l];
			cout << fibers_length[l];

			nodes_on_fibers[l] = (long int *)malloc(sizeof(long int) * fibers_length[l]);
			for (i = 0; i < fibers_length[l]; i++) {
				fibers_vtkfile >> nodes_on_fibers[l][i];
				cout << " " << nodes_on_fibers[l][i];
			}
			cout << endl;
		}

		long int count = 0;

		for (i = 0; i < Nfiber; i++)
			for (j = 0; j < fibers_length[i]; j++) {
				node_index = nodes_on_fibers[i][j];
				nodes[node_index]->on_fibers = count;
				count++;
			}
	}

	fibers_vtkfile.close();

	cout << "Initialized fibers" << endl;
	// Done //////////////////////////////////////////////////////////////////

	//Read late xyz and update coordinates
	ifstream latevtkfile("late169.gel.vtk");
	if (!latevtkfile.is_open())
	{
			cout << "Error opening late vtk file" << endl;
			exit(1);
	}
	char latebuffer[256];
	latevtkfile.getline(latebuffer, 100);
	latevtkfile.getline(latebuffer, 100);
	latevtkfile.getline(latebuffer, 100);
	latevtkfile.getline(latebuffer, 100);
	latevtkfile >> latebuffer;
	latevtkfile >> N_nodes;
	latevtkfile >> latebuffer;
	for (i = 0; i < N_nodes; i++){
	latevtkfile >> nodes[i]->position[0]>> nodes[i]->position[1]
			>>nodes[i]->position[2];
	}
		////////////////////////////////////////////////////////////////////////////////

	latevtkfile.close();
	cout<<"finish reading xyz coordinate from late stage"<<endl;
	return 0;
}


void CalculateVolumes() {
	for (long int element_index = 0; element_index < N_elements; element_index++) {
		elements[element_index]->update_volume();
		elements[element_index]->phi = Vol0_phi0 / elements[element_index]->volume;
	}
}

void CalculatePressure(double t0) {
	double phi_tmp;
	double vSP_tmp;

	if (t0 >= TimeLightON)
		LightSpace = ILIGHT;
	else
		LightSpace = 0.;

	for (long int element_index = 0; element_index < N_elements; element_index++) {
		phi_tmp = elements[element_index]->phi;
		vSP_tmp = elements[element_index]->vSP;
		elements[element_index]->vSP = vSP_tmp + (LightSpace * (elements[element_index]->vSPTOT - vSP_tmp) - Kdark * vSP_tmp) * h;
		elements[element_index]->P0 = -(phi_tmp + log(1.0 - phi_tmp) + (chi0 + chi1 * phi_tmp) * phi_tmp * phi_tmp)
			+ c0 *Vol0 * phi_tmp / (2.0 * phi0) - kappa * phi_tmp * vSP_tmp;
	}
}


void CalculateMobilities() {
	double phi_aver;

	for (long int node_index = 0; node_index < N_nodes; node_index++) {
		phi_aver = nodes[node_index]->compute_phi();
		nodes[node_index]->mobility = L0 * (1.0 - phi_aver) / sqrt(phi_aver) * 8.0;
		if (MobilityonTemperatureON == 1)
			nodes[node_index]->mobility = nodes[node_index]->mobility * (temperature + 273.) / (22. + 273.);
	}

}


//Ya's changes: 1. if nodes on bump, velocity = 0
//2. if nodes on the side surface, nodes move laterally
void isBump(){

	for (int node_index = 0; node_index < N_nodes; node_index++){
		int px=nodes[node_index]->position[0]/unit_space;
		if(px>Lx)
			px = px +1;

		int py=nodes[node_index]->position[1]/unit_space;
		if(py>nodes[node_index]->position[1]/unit_space-0.2)
			py = py +1;
		int pz=nodes[node_index]->position[2]/unit_space;
		if(pz>nodes[node_index]->position[2]/unit_space-0.2)
			pz = pz +1;

		double bumpedge = sqrt(bump_radius*bump_radius-(px-0.5*Lx)*(px-0.5*Lx))
			   -(bump_radius-bump_height)-1;
		//cout<<px<<'\t'<<py<<'\t'<<pz<<'\t'<<bumpedge<<'\t'<<Ly<<endl;
		if(py == 0 || py == Ly){ // check if is edge surface.unit_space
			nodes[node_index]->Pos_Node = 1;
				 //   cout<<"edge node"<<endl;
				 }
		else if(pz < bumpedge){
	    	nodes[node_index]->Pos_Node = 0;
		}
		else{
			nodes[node_index]->Pos_Node = 2;
		//    cout<<"body node"<<endl;
		}
	}
}

void CalculateVelocities() {
	long int node_index;
	long int element_index, node_index_down, node_index_up;
	long int i, j, k;
	double   Mob0, lsprx, lspry, lsprz, lengthtopspr, Fabs, deltaTheta;
	double   rxD, rxU, ryD, ryU, rzD, rzU, lD, lU, dotprod, cosT, sinT,
		FabsBendU, FabsBendD, thetaT, FbendFactor;
	double fbUX, fbUY, fbUZ, fbDX, fbDY, fbDZ;

	rxD = 0.; rxU = 0.; ryD = 0.; ryU = 0.; FbendFactor = 0.;
	rzD = 0.; rzU = 0.; lD = 0.; lU = 0.; dotprod = 0.; cosT = 0.;
	FabsBendU = 0.;
	FabsBendD = 0.; thetaT = 3.1415;

	for (node_index = 0; node_index < N_nodes; node_index++) {
		nodes[node_index]->velocity[0] = 0.0;
		nodes[node_index]->velocity[1] = 0.0;
		nodes[node_index]->velocity[2] = 0.0;
	}

	for (node_index = 0; node_index < N_nodes; node_index++) {
		nodes[node_index]->FSpring[0] = 0.0;
		nodes[node_index]->FSpring[1] = 0.0;
		nodes[node_index]->FSpring[2] = 0.0;
		nodes[node_index]->FPressure[0] = 0.0;
		nodes[node_index]->FPressure[1] = 0.0;
		nodes[node_index]->FPressure[2] = 0.0;
		nodes[node_index]->FCilia[0] = 0.0;
		nodes[node_index]->FCilia[1] = 0.0;
		nodes[node_index]->FCilia[2] = 0.0;
		nodes[node_index]->FCiliaBend[0] = 0.0;
		nodes[node_index]->FCiliaBend[1] = 0.0;
		nodes[node_index]->FCiliaBend[2] = 0.0;
	}

	for (element_index = 0; element_index < N_elements; element_index++) {
		elements[element_index]->update_nodes_FSpring();
		elements[element_index]->update_nodes_FPressure();
	}

	// Add Hookean Springs: notations are on
	for (i = 0; i < Nfiber; i++)
		for (j = 0; j < fibers_length[i] - 1; j++) {
			node_index = nodes_on_fibers[i][j];
			node_index_up = nodes_on_fibers[i][j + 1];
			lsprx = nodes[node_index]->position[0] - nodes[node_index_up]->position[0];
			lspry = nodes[node_index]->position[1] - nodes[node_index_up]->position[1];
			lsprz = nodes[node_index]->position[2] - nodes[node_index_up]->position[2];
			lengthtopspr = sqrt(lsprx * lsprx + lspry * lspry + lsprz * lsprz);
			Fabs = Bcilia * (-1. + fiber_equi_length / lengthtopspr);
			nodes[node_index]->FCilia[0] = nodes[node_index]->FCilia[0] + Fabs * lsprx;
			nodes[node_index_up]->FCilia[0] = nodes[node_index_up]->FCilia[0] - Fabs * lsprx;
			nodes[node_index]->FCilia[1] = nodes[node_index]->FCilia[1] + Fabs * lspry;
			nodes[node_index_up]->FCilia[1] = nodes[node_index_up]->FCilia[1] - Fabs * lspry;
			nodes[node_index]->FCilia[2] = nodes[node_index]->FCilia[2] + Fabs * lsprz;
			nodes[node_index_up]->FCilia[2] = nodes[node_index_up]->FCilia[2] - Fabs * lsprz;
		}

	// for angular potential
	for (i = 0; i < Nfiber; i++)
		for (j = 0; j < fibers_length[i] - 2; j++) {
			node_index = nodes_on_fibers[i][j + 1];
			node_index_down = nodes_on_fibers[i][j];
			node_index_up = nodes_on_fibers[i][j + 2];
			rxD = nodes[node_index_down]->position[0] - nodes[node_index]->position[0];
			ryD = nodes[node_index_down]->position[1] - nodes[node_index]->position[1];
			rzD = nodes[node_index_down]->position[2] - nodes[node_index]->position[2];
			rxU = nodes[node_index_up]->position[0] - nodes[node_index]->position[0];
			ryU = nodes[node_index_up]->position[1] - nodes[node_index]->position[1];
			rzU = nodes[node_index_up]->position[2] - nodes[node_index]->position[2];
			lD = sqrt(rxD * rxD + ryD * ryD + rzD * rzD);
			lU = sqrt(rxU * rxU + ryU * ryU + rzU * rzU);
			dotprod = rxD * rxU + ryD * ryU + rzD * rzU;
			cosT = dotprod / (lD * lU);
			thetaT = acos(cosT);
			sinT = sin(thetaT);
			deltaTheta = thetaT - 3.14159265358979323846264;

			if (fabs(deltaTheta) < 0.0000000001)
				FbendFactor = Kbend * (-1. - deltaTheta * deltaTheta / 6.);
			else
				FbendFactor = Kbend * deltaTheta / sinT;

			if (UseCosHarmonicBend == 1)
				FbendFactor = -1. * Kbend * (cosT + (double)1.);

			fbUX = (cosT * rxU / lU - rxD / lD) * FbendFactor / lU;
			fbUY = (cosT * ryU / lU - ryD / lD) * FbendFactor / lU;
			fbUZ = (cosT * rzU / lU - rzD / lD) * FbendFactor / lU;
			fbDX = (cosT * rxD / lD - rxU / lU) * FbendFactor / lD;
			fbDY = (cosT * ryD / lD - ryU / lU) * FbendFactor / lD;
			fbDZ = (cosT * rzD / lD - rzU / lU) * FbendFactor / lD;

			nodes[node_index]->FCiliaBend[0] = nodes[node_index]->FCiliaBend[0] + (fbUX + fbDX);
			nodes[node_index]->FCiliaBend[1] = nodes[node_index]->FCiliaBend[1] + (fbUY + fbDY);
			nodes[node_index]->FCiliaBend[2] = nodes[node_index]->FCiliaBend[2] + (fbUZ + fbDZ);

			nodes[node_index_up]->FCiliaBend[0] = nodes[node_index_up]->FCiliaBend[0] - fbUX;
			nodes[node_index_up]->FCiliaBend[1] = nodes[node_index_up]->FCiliaBend[1] - fbUY;
			nodes[node_index_up]->FCiliaBend[2] = nodes[node_index_up]->FCiliaBend[2] - fbUZ;

			nodes[node_index_down]->FCiliaBend[0] = nodes[node_index_down]->FCiliaBend[0] - fbDX;
			nodes[node_index_down]->FCiliaBend[1] = nodes[node_index_down]->FCiliaBend[1] - fbDY;
			nodes[node_index_down]->FCiliaBend[2] = nodes[node_index_down]->FCiliaBend[2] - fbDZ;
		}

	// update velocities
	//Ya's change: depend on location
	for (node_index = 0; node_index < N_nodes; node_index++) {
		if(nodes[node_index]->Pos_Node==0){
			nodes[node_index]->velocity[0] = 0;
			nodes[node_index]->velocity[1] = 0;
			nodes[node_index]->velocity[2] = 0;
		}
		else if(nodes[node_index]->Pos_Node==1 || nodes[node_index]->Pos_Node==2){
		nodes[node_index]->velocity[0] = nodes[node_index]->mobility *
			(nodes[node_index]->FSpring[0] + nodes[node_index]->FPressure[0] + nodes[node_index]->FCilia[0] + nodes[node_index]->FCiliaBend[0]);
		nodes[node_index]->velocity[1] = 0;
		nodes[node_index]->velocity[2] = nodes[node_index]->mobility *
			(nodes[node_index]->FSpring[2] + nodes[node_index]->FPressure[2] + nodes[node_index]->FCilia[2] + nodes[node_index]->FCiliaBend[2]);

		}
		else{
		nodes[node_index]->velocity[0] = nodes[node_index]->mobility *
			(nodes[node_index]->FSpring[0] + nodes[node_index]->FPressure[0] + nodes[node_index]->FCilia[0] + nodes[node_index]->FCiliaBend[0]);
		nodes[node_index]->velocity[1] = nodes[node_index]->mobility *
			(nodes[node_index]->FSpring[1] + nodes[node_index]->FPressure[1] + nodes[node_index]->FCilia[1] + nodes[node_index]->FCiliaBend[1]);
		nodes[node_index]->velocity[2] = nodes[node_index]->mobility *
			(nodes[node_index]->FSpring[2] + nodes[node_index]->FPressure[2] + nodes[node_index]->FCilia[2] + nodes[node_index]->FCiliaBend[2]);
		}
	}
}


// update coordinate: isBump = 0:no update; =1: update x,z; =2 update x,y,z
void Update_XYZ() {
	long int node_index;
	for (node_index = 0; node_index < N_nodes; node_index++) {
		nodes[node_index]->position[0] = nodes[node_index]->position[0] + nodes[node_index]->velocity[0] * h;
		nodes[node_index]->position[1] = nodes[node_index]->position[1] + nodes[node_index]->velocity[1] * h;
		nodes[node_index]->position[2] = nodes[node_index]->position[2] + nodes[node_index]->velocity[2] * h;
	}
}

void Element::update_volume() {
	double tmp, X[8], Y[8], Z[8];
	int    i;
	long int node_index;

	for (i = 0; i < 8; i++)
	{
		node_index = nodes_index[i];
		X[i] = nodes_array[node_index]->position[0];
		Y[i] = nodes_array[node_index]->position[1];
		Z[i] = nodes_array[node_index]->position[2];
	}

	tmp = X[2] * Y[1] * Z[0] + X[3] * Y[1] * Z[0] - X[4] * Y[1] * Z[0]
		- X[5] * Y[1] * Z[0] - X[1] * Y[2] * Z[0] + X[3] * Y[2] * Z[0]
		- X[1] * Y[3] * Z[0] - X[2] * Y[3] * Z[0] + X[4] * Y[3] * Z[0]
		+ X[7] * Y[3] * Z[0] + X[1] * Y[4] * Z[0] - X[3] * Y[4] * Z[0]
		+ X[5] * Y[4] * Z[0] - X[7] * Y[4] * Z[0] + X[1] * Y[5] * Z[0]
		- X[4] * Y[5] * Z[0] - X[3] * Y[7] * Z[0] + X[4] * Y[7] * Z[0]
		- X[2] * Y[0] * Z[1] - X[3] * Y[0] * Z[1] + X[4] * Y[0] * Z[1]
		+ X[5] * Y[0] * Z[1] + X[0] * Y[2] * Z[1] + X[3] * Y[2] * Z[1]
		- X[5] * Y[2] * Z[1] - X[6] * Y[2] * Z[1] + X[0] * Y[3] * Z[1]
		- X[2] * Y[3] * Z[1] - X[0] * Y[4] * Z[1] + X[5] * Y[4] * Z[1]
		- X[0] * Y[5] * Z[1] + X[2] * Y[5] * Z[1] - X[4] * Y[5] * Z[1]
		+ X[6] * Y[5] * Z[1] + X[2] * Y[6] * Z[1] - X[5] * Y[6] * Z[1]
		+ X[1] * Y[0] * Z[2] - X[3] * Y[0] * Z[2] - X[0] * Y[1] * Z[2]
		- X[3] * Y[1] * Z[2] + X[5] * Y[1] * Z[2] + X[6] * Y[1] * Z[2]
		+ X[0] * Y[3] * Z[2] + X[1] * Y[3] * Z[2] - X[6] * Y[3] * Z[2]
		- X[7] * Y[3] * Z[2] - X[1] * Y[5] * Z[2] + X[6] * Y[5] * Z[2]
		- X[1] * Y[6] * Z[2] + X[3] * Y[6] * Z[2] - X[5] * Y[6] * Z[2]
		+ X[7] * Y[6] * Z[2] + X[3] * Y[7] * Z[2] - X[6] * Y[7] * Z[2]
		+ X[1] * Y[0] * Z[3] + X[2] * Y[0] * Z[3] - X[4] * Y[0] * Z[3]
		- X[7] * Y[0] * Z[3] - X[0] * Y[1] * Z[3] + X[2] * Y[1] * Z[3]
		- X[0] * Y[2] * Z[3] - X[1] * Y[2] * Z[3] + X[6] * Y[2] * Z[3]
		+ X[7] * Y[2] * Z[3] + X[0] * Y[4] * Z[3] - X[7] * Y[4] * Z[3]
		- X[2] * Y[6] * Z[3] + X[7] * Y[6] * Z[3] + X[0] * Y[7] * Z[3]
		- X[2] * Y[7] * Z[3] + X[4] * Y[7] * Z[3] - X[6] * Y[7] * Z[3]
		- X[1] * Y[0] * Z[4] + X[3] * Y[0] * Z[4] - X[5] * Y[0] * Z[4]
		+ X[7] * Y[0] * Z[4] + X[0] * Y[1] * Z[4] - X[5] * Y[1] * Z[4]
		- X[0] * Y[3] * Z[4] + X[7] * Y[3] * Z[4] + X[0] * Y[5] * Z[4]
		+ X[1] * Y[5] * Z[4] - X[6] * Y[5] * Z[4] - X[7] * Y[5] * Z[4]
		+ X[5] * Y[6] * Z[4] - X[7] * Y[6] * Z[4] - X[0] * Y[7] * Z[4]
		- X[3] * Y[7] * Z[4] + X[5] * Y[7] * Z[4] + X[6] * Y[7] * Z[4]
		- X[1] * Y[0] * Z[5] + X[4] * Y[0] * Z[5] + X[0] * Y[1] * Z[5]
		- X[2] * Y[1] * Z[5] + X[4] * Y[1] * Z[5] - X[6] * Y[1] * Z[5]
		+ X[1] * Y[2] * Z[5] - X[6] * Y[2] * Z[5] - X[0] * Y[4] * Z[5]
		- X[1] * Y[4] * Z[5] + X[6] * Y[4] * Z[5] + X[7] * Y[4] * Z[5]
		+ X[1] * Y[6] * Z[5] + X[2] * Y[6] * Z[5] - X[4] * Y[6] * Z[5]
		- X[7] * Y[6] * Z[5] - X[4] * Y[7] * Z[5] + X[6] * Y[7] * Z[5]
		- X[2] * Y[1] * Z[6] + X[5] * Y[1] * Z[6] + X[1] * Y[2] * Z[6]
		- X[3] * Y[2] * Z[6] + X[5] * Y[2] * Z[6] - X[7] * Y[2] * Z[6]
		+ X[2] * Y[3] * Z[6] - X[7] * Y[3] * Z[6] - X[5] * Y[4] * Z[6]
		+ X[7] * Y[4] * Z[6] - X[1] * Y[5] * Z[6] - X[2] * Y[5] * Z[6]
		+ X[4] * Y[5] * Z[6] + X[7] * Y[5] * Z[6] + X[2] * Y[7] * Z[6]
		+ X[3] * Y[7] * Z[6] - X[4] * Y[7] * Z[6] - X[5] * Y[7] * Z[6]
		+ X[3] * Y[0] * Z[7] - X[4] * Y[0] * Z[7] - X[3] * Y[2] * Z[7]
		+ X[6] * Y[2] * Z[7] - X[0] * Y[3] * Z[7] + X[2] * Y[3] * Z[7]
		- X[4] * Y[3] * Z[7] + X[6] * Y[3] * Z[7] + X[0] * Y[4] * Z[7]
		+ X[3] * Y[4] * Z[7] - X[5] * Y[4] * Z[7] - X[6] * Y[4] * Z[7]
		+ X[4] * Y[5] * Z[7] - X[6] * Y[5] * Z[7] - X[2] * Y[6] * Z[7]
		- X[3] * Y[6] * Z[7] + X[4] * Y[6] * Z[7] + X[5] * Y[6] * Z[7];

	tmp = tmp / 12.;
	volume = tmp;
}

void Element::update_nodes_FSpring() {

	double x13 = nodes_array[nodes_index[2]]->position[0] - nodes_array[nodes_index[0]]->position[0];
	double x24 = nodes_array[nodes_index[3]]->position[0] - nodes_array[nodes_index[1]]->position[0];
	double x57 = nodes_array[nodes_index[6]]->position[0] - nodes_array[nodes_index[4]]->position[0];
	double x68 = nodes_array[nodes_index[7]]->position[0] - nodes_array[nodes_index[5]]->position[0];
	double x18 = nodes_array[nodes_index[7]]->position[0] - nodes_array[nodes_index[0]]->position[0];
	double x45 = nodes_array[nodes_index[4]]->position[0] - nodes_array[nodes_index[3]]->position[0];
	double x27 = nodes_array[nodes_index[6]]->position[0] - nodes_array[nodes_index[1]]->position[0];
	double x36 = nodes_array[nodes_index[5]]->position[0] - nodes_array[nodes_index[2]]->position[0];
	double x16 = nodes_array[nodes_index[5]]->position[0] - nodes_array[nodes_index[0]]->position[0];
	double x25 = nodes_array[nodes_index[4]]->position[0] - nodes_array[nodes_index[1]]->position[0];
	double x47 = nodes_array[nodes_index[6]]->position[0] - nodes_array[nodes_index[3]]->position[0];
	double x38 = nodes_array[nodes_index[7]]->position[0] - nodes_array[nodes_index[2]]->position[0];
	double x17 = nodes_array[nodes_index[6]]->position[0] - nodes_array[nodes_index[0]]->position[0];
	double x35 = nodes_array[nodes_index[4]]->position[0] - nodes_array[nodes_index[2]]->position[0];
	double x46 = nodes_array[nodes_index[5]]->position[0] - nodes_array[nodes_index[3]]->position[0];
	double x28 = nodes_array[nodes_index[7]]->position[0] - nodes_array[nodes_index[1]]->position[0];

	double y13 = nodes_array[nodes_index[2]]->position[1] - nodes_array[nodes_index[0]]->position[1];
	double y24 = nodes_array[nodes_index[3]]->position[1] - nodes_array[nodes_index[1]]->position[1];
	double y57 = nodes_array[nodes_index[6]]->position[1] - nodes_array[nodes_index[4]]->position[1];
	double y68 = nodes_array[nodes_index[7]]->position[1] - nodes_array[nodes_index[5]]->position[1];
	double y18 = nodes_array[nodes_index[7]]->position[1] - nodes_array[nodes_index[0]]->position[1];
	double y45 = nodes_array[nodes_index[4]]->position[1] - nodes_array[nodes_index[3]]->position[1];
	double y27 = nodes_array[nodes_index[6]]->position[1] - nodes_array[nodes_index[1]]->position[1];
	double y36 = nodes_array[nodes_index[5]]->position[1] - nodes_array[nodes_index[2]]->position[1];
	double y16 = nodes_array[nodes_index[5]]->position[1] - nodes_array[nodes_index[0]]->position[1];
	double y25 = nodes_array[nodes_index[4]]->position[1] - nodes_array[nodes_index[1]]->position[1];
	double y47 = nodes_array[nodes_index[6]]->position[1] - nodes_array[nodes_index[3]]->position[1];
	double y38 = nodes_array[nodes_index[7]]->position[1] - nodes_array[nodes_index[2]]->position[1];
	double y17 = nodes_array[nodes_index[6]]->position[1] - nodes_array[nodes_index[0]]->position[1];
	double y35 = nodes_array[nodes_index[4]]->position[1] - nodes_array[nodes_index[2]]->position[1];
	double y46 = nodes_array[nodes_index[5]]->position[1] - nodes_array[nodes_index[3]]->position[1];
	double y28 = nodes_array[nodes_index[7]]->position[1] - nodes_array[nodes_index[1]]->position[1];

	double z13 = nodes_array[nodes_index[2]]->position[2] - nodes_array[nodes_index[0]]->position[2];
	double z24 = nodes_array[nodes_index[3]]->position[2] - nodes_array[nodes_index[1]]->position[2];
	double z57 = nodes_array[nodes_index[6]]->position[2] - nodes_array[nodes_index[4]]->position[2];
	double z68 = nodes_array[nodes_index[7]]->position[2] - nodes_array[nodes_index[5]]->position[2];
	double z18 = nodes_array[nodes_index[7]]->position[2] - nodes_array[nodes_index[0]]->position[2];
	double z45 = nodes_array[nodes_index[4]]->position[2] - nodes_array[nodes_index[3]]->position[2];
	double z27 = nodes_array[nodes_index[6]]->position[2] - nodes_array[nodes_index[1]]->position[2];
	double z36 = nodes_array[nodes_index[5]]->position[2] - nodes_array[nodes_index[2]]->position[2];
	double z16 = nodes_array[nodes_index[5]]->position[2] - nodes_array[nodes_index[0]]->position[2];
	double z25 = nodes_array[nodes_index[4]]->position[2] - nodes_array[nodes_index[1]]->position[2];
	double z47 = nodes_array[nodes_index[6]]->position[2] - nodes_array[nodes_index[3]]->position[2];
	double z38 = nodes_array[nodes_index[7]]->position[2] - nodes_array[nodes_index[2]]->position[2];
	double z17 = nodes_array[nodes_index[6]]->position[2] - nodes_array[nodes_index[0]]->position[2];
	double z35 = nodes_array[nodes_index[4]]->position[2] - nodes_array[nodes_index[2]]->position[2];
	double z46 = nodes_array[nodes_index[5]]->position[2] - nodes_array[nodes_index[3]]->position[2];
	double z28 = nodes_array[nodes_index[7]]->position[2] - nodes_array[nodes_index[1]]->position[2];

	//13
	nodes_array[nodes_index[0]]->FSpring[0] = nodes_array[nodes_index[0]]->FSpring[0]
		+ c06*x13;
	nodes_array[nodes_index[2]]->FSpring[0] = nodes_array[nodes_index[2]]->FSpring[0]
		- c06*x13;
	//24
	nodes_array[nodes_index[1]]->FSpring[0] = nodes_array[nodes_index[1]]->FSpring[0]
		+ c06*x24;
	nodes_array[nodes_index[3]]->FSpring[0] = nodes_array[nodes_index[3]]->FSpring[0]
		- c06*x24;
	//57
	nodes_array[nodes_index[4]]->FSpring[0] = nodes_array[nodes_index[4]]->FSpring[0]
		+ c06*x57;
	nodes_array[nodes_index[6]]->FSpring[0] = nodes_array[nodes_index[6]]->FSpring[0]
		- c06*x57;
	//68
	nodes_array[nodes_index[5]]->FSpring[0] = nodes_array[nodes_index[5]]->FSpring[0]
		+ c06*x68;
	nodes_array[nodes_index[7]]->FSpring[0] = nodes_array[nodes_index[7]]->FSpring[0]
		- c06*x68;
	//18
	nodes_array[nodes_index[0]]->FSpring[0] = nodes_array[nodes_index[0]]->FSpring[0]
		+ c06*x18;
	nodes_array[nodes_index[7]]->FSpring[0] = nodes_array[nodes_index[7]]->FSpring[0]
		- c06*x18;
	//45
	nodes_array[nodes_index[3]]->FSpring[0] = nodes_array[nodes_index[3]]->FSpring[0]
		+ c06*x45;
	nodes_array[nodes_index[4]]->FSpring[0] = nodes_array[nodes_index[4]]->FSpring[0]
		- c06*x45;
	//27
	nodes_array[nodes_index[1]]->FSpring[0] = nodes_array[nodes_index[1]]->FSpring[0]
		+ c06*x27;
	nodes_array[nodes_index[6]]->FSpring[0] = nodes_array[nodes_index[6]]->FSpring[0]
		- c06*x27;
	//36
	nodes_array[nodes_index[2]]->FSpring[0] = nodes_array[nodes_index[2]]->FSpring[0]
		+ c06*x36;
	nodes_array[nodes_index[5]]->FSpring[0] = nodes_array[nodes_index[5]]->FSpring[0]
		- c06*x36;
	//16
	nodes_array[nodes_index[0]]->FSpring[0] = nodes_array[nodes_index[0]]->FSpring[0]
		+ c06*x16;
	nodes_array[nodes_index[5]]->FSpring[0] = nodes_array[nodes_index[5]]->FSpring[0]
		- c06*x16;
	//25
	nodes_array[nodes_index[1]]->FSpring[0] = nodes_array[nodes_index[1]]->FSpring[0]
		+ c06*x25;
	nodes_array[nodes_index[4]]->FSpring[0] = nodes_array[nodes_index[4]]->FSpring[0]
		- c06*x25;
	//47
	nodes_array[nodes_index[3]]->FSpring[0] = nodes_array[nodes_index[3]]->FSpring[0]
		+ c06*x47;
	nodes_array[nodes_index[6]]->FSpring[0] = nodes_array[nodes_index[6]]->FSpring[0]
		- c06*x47;
	//38
	nodes_array[nodes_index[2]]->FSpring[0] = nodes_array[nodes_index[2]]->FSpring[0]
		+ c06*x38;
	nodes_array[nodes_index[7]]->FSpring[0] = nodes_array[nodes_index[7]]->FSpring[0]
		- c06*x38;
	//17
	nodes_array[nodes_index[0]]->FSpring[0] = nodes_array[nodes_index[0]]->FSpring[0]
		+ c06*x17;
	nodes_array[nodes_index[6]]->FSpring[0] = nodes_array[nodes_index[6]]->FSpring[0]
		- c06*x17;
	//35
	nodes_array[nodes_index[2]]->FSpring[0] = nodes_array[nodes_index[2]]->FSpring[0]
		+ c06*x35;
	nodes_array[nodes_index[4]]->FSpring[0] = nodes_array[nodes_index[4]]->FSpring[0]
		- c06*x35;
	//46
	nodes_array[nodes_index[3]]->FSpring[0] = nodes_array[nodes_index[3]]->FSpring[0]
		+ c06*x46;
	nodes_array[nodes_index[5]]->FSpring[0] = nodes_array[nodes_index[5]]->FSpring[0]
		- c06*x46;
	//28
	nodes_array[nodes_index[1]]->FSpring[0] = nodes_array[nodes_index[1]]->FSpring[0]
		+ c06*x28;
	nodes_array[nodes_index[7]]->FSpring[0] = nodes_array[nodes_index[7]]->FSpring[0]
		- c06*x28;

	//13
	nodes_array[nodes_index[0]]->FSpring[1] = nodes_array[nodes_index[0]]->FSpring[1]
		+ c06*y13;
	nodes_array[nodes_index[2]]->FSpring[1] = nodes_array[nodes_index[2]]->FSpring[1]
		- c06*y13;
	//24
	nodes_array[nodes_index[1]]->FSpring[1] = nodes_array[nodes_index[1]]->FSpring[1]
		+ c06*y24;
	nodes_array[nodes_index[3]]->FSpring[1] = nodes_array[nodes_index[3]]->FSpring[1]
		- c06*y24;
	//57
	nodes_array[nodes_index[4]]->FSpring[1] = nodes_array[nodes_index[4]]->FSpring[1]
		+ c06*y57;
	nodes_array[nodes_index[6]]->FSpring[1] = nodes_array[nodes_index[6]]->FSpring[1]
		- c06*y57;
	//68
	nodes_array[nodes_index[5]]->FSpring[1] = nodes_array[nodes_index[5]]->FSpring[1]
		+ c06*y68;
	nodes_array[nodes_index[7]]->FSpring[1] = nodes_array[nodes_index[7]]->FSpring[1]
		- c06*y68;
	//18
	nodes_array[nodes_index[0]]->FSpring[1] = nodes_array[nodes_index[0]]->FSpring[1]
		+ c06*y18;
	nodes_array[nodes_index[7]]->FSpring[1] = nodes_array[nodes_index[7]]->FSpring[1]
		- c06*y18;
	//45
	nodes_array[nodes_index[3]]->FSpring[1] = nodes_array[nodes_index[3]]->FSpring[1]
		+ c06*y45;
	nodes_array[nodes_index[4]]->FSpring[1] = nodes_array[nodes_index[4]]->FSpring[1]
		- c06*y45;
	//27
	nodes_array[nodes_index[1]]->FSpring[1] = nodes_array[nodes_index[1]]->FSpring[1]
		+ c06*y27;
	nodes_array[nodes_index[6]]->FSpring[1] = nodes_array[nodes_index[6]]->FSpring[1]
		- c06*y27;
	//36
	nodes_array[nodes_index[2]]->FSpring[1] = nodes_array[nodes_index[2]]->FSpring[1]
		+ c06*y36;
	nodes_array[nodes_index[5]]->FSpring[1] = nodes_array[nodes_index[5]]->FSpring[1]
		- c06*y36;
	//16
	nodes_array[nodes_index[0]]->FSpring[1] = nodes_array[nodes_index[0]]->FSpring[1]
		+ c06*y16;
	nodes_array[nodes_index[5]]->FSpring[1] = nodes_array[nodes_index[5]]->FSpring[1]
		- c06*y16;
	//25
	nodes_array[nodes_index[1]]->FSpring[1] = nodes_array[nodes_index[1]]->FSpring[1]
		+ c06*y25;
	nodes_array[nodes_index[4]]->FSpring[1] = nodes_array[nodes_index[4]]->FSpring[1]
		- c06*y25;
	//47
	nodes_array[nodes_index[3]]->FSpring[1] = nodes_array[nodes_index[3]]->FSpring[1]
		+ c06*y47;
	nodes_array[nodes_index[6]]->FSpring[1] = nodes_array[nodes_index[6]]->FSpring[1]
		- c06*y47;
	//38
	nodes_array[nodes_index[2]]->FSpring[1] = nodes_array[nodes_index[2]]->FSpring[1]
		+ c06*y38;
	nodes_array[nodes_index[7]]->FSpring[1] = nodes_array[nodes_index[7]]->FSpring[1]
		- c06*y38;
	//17
	nodes_array[nodes_index[0]]->FSpring[1] = nodes_array[nodes_index[0]]->FSpring[1]
		+ c06*y17;
	nodes_array[nodes_index[6]]->FSpring[1] = nodes_array[nodes_index[6]]->FSpring[1]
		- c06*y17;
	//35
	nodes_array[nodes_index[2]]->FSpring[1] = nodes_array[nodes_index[2]]->FSpring[1]
		+ c06*y35;
	nodes_array[nodes_index[4]]->FSpring[1] = nodes_array[nodes_index[4]]->FSpring[1]
		- c06*y35;
	//46
	nodes_array[nodes_index[3]]->FSpring[1] = nodes_array[nodes_index[3]]->FSpring[1]
		+ c06*y46;
	nodes_array[nodes_index[5]]->FSpring[1] = nodes_array[nodes_index[5]]->FSpring[1]
		- c06*y46;
	//28
	nodes_array[nodes_index[1]]->FSpring[1] = nodes_array[nodes_index[1]]->FSpring[1]
		+ c06*y28;
	nodes_array[nodes_index[7]]->FSpring[1] = nodes_array[nodes_index[7]]->FSpring[1]
		- c06*y28;

	//13
	nodes_array[nodes_index[0]]->FSpring[2] = nodes_array[nodes_index[0]]->FSpring[2]
		+ c06*z13;
	nodes_array[nodes_index[2]]->FSpring[2] = nodes_array[nodes_index[2]]->FSpring[2]
		- c06*z13;
	//24
	nodes_array[nodes_index[1]]->FSpring[2] = nodes_array[nodes_index[1]]->FSpring[2]
		+ c06*z24;
	nodes_array[nodes_index[3]]->FSpring[2] = nodes_array[nodes_index[3]]->FSpring[2]
		- c06*z24;
	//57
	nodes_array[nodes_index[4]]->FSpring[2] = nodes_array[nodes_index[4]]->FSpring[2]
		+ c06*z57;
	nodes_array[nodes_index[6]]->FSpring[2] = nodes_array[nodes_index[6]]->FSpring[2]
		- c06*z57;
	//68
	nodes_array[nodes_index[5]]->FSpring[2] = nodes_array[nodes_index[5]]->FSpring[2]
		+ c06*z68;
	nodes_array[nodes_index[7]]->FSpring[2] = nodes_array[nodes_index[7]]->FSpring[2]
		- c06*z68;
	//18
	nodes_array[nodes_index[0]]->FSpring[2] = nodes_array[nodes_index[0]]->FSpring[2]
		+ c06*z18;
	nodes_array[nodes_index[7]]->FSpring[2] = nodes_array[nodes_index[7]]->FSpring[2]
		- c06*z18;
	//45
	nodes_array[nodes_index[3]]->FSpring[2] = nodes_array[nodes_index[3]]->FSpring[2]
		+ c06*z45;
	nodes_array[nodes_index[4]]->FSpring[2] = nodes_array[nodes_index[4]]->FSpring[2]
		- c06*z45;
	//27
	nodes_array[nodes_index[1]]->FSpring[2] = nodes_array[nodes_index[1]]->FSpring[2]
		+ c06*z27;
	nodes_array[nodes_index[6]]->FSpring[2] = nodes_array[nodes_index[6]]->FSpring[2]
		- c06*z27;
	//36
	nodes_array[nodes_index[2]]->FSpring[2] = nodes_array[nodes_index[2]]->FSpring[2]
		+ c06*z36;
	nodes_array[nodes_index[5]]->FSpring[2] = nodes_array[nodes_index[5]]->FSpring[2]
		- c06*z36;
	//16
	nodes_array[nodes_index[0]]->FSpring[2] = nodes_array[nodes_index[0]]->FSpring[2]
		+ c06*z16;
	nodes_array[nodes_index[5]]->FSpring[2] = nodes_array[nodes_index[5]]->FSpring[2]
		- c06*z16;
	//25
	nodes_array[nodes_index[1]]->FSpring[2] = nodes_array[nodes_index[1]]->FSpring[2]
		+ c06*z25;
	nodes_array[nodes_index[4]]->FSpring[2] = nodes_array[nodes_index[4]]->FSpring[2]
		- c06*z25;
	//47
	nodes_array[nodes_index[3]]->FSpring[2] = nodes_array[nodes_index[3]]->FSpring[2]
		+ c06*z47;
	nodes_array[nodes_index[6]]->FSpring[2] = nodes_array[nodes_index[6]]->FSpring[2]
		- c06*z47;
	//38
	nodes_array[nodes_index[2]]->FSpring[2] = nodes_array[nodes_index[2]]->FSpring[2]
		+ c06*z38;
	nodes_array[nodes_index[7]]->FSpring[2] = nodes_array[nodes_index[7]]->FSpring[2]
		- c06*z38;
	//17
	nodes_array[nodes_index[0]]->FSpring[2] = nodes_array[nodes_index[0]]->FSpring[2]
		+ c06*z17;
	nodes_array[nodes_index[6]]->FSpring[2] = nodes_array[nodes_index[6]]->FSpring[2]
		- c06*z17;
	//35
	nodes_array[nodes_index[2]]->FSpring[2] = nodes_array[nodes_index[2]]->FSpring[2]
		+ c06*z35;
	nodes_array[nodes_index[4]]->FSpring[2] = nodes_array[nodes_index[4]]->FSpring[2]
		- c06*z35;
	//46
	nodes_array[nodes_index[3]]->FSpring[2] = nodes_array[nodes_index[3]]->FSpring[2]
		+ c06*z46;
	nodes_array[nodes_index[5]]->FSpring[2] = nodes_array[nodes_index[5]]->FSpring[2]
		- c06*z46;
	//28
	nodes_array[nodes_index[1]]->FSpring[2] = nodes_array[nodes_index[1]]->FSpring[2]
		+ c06*z28;
	nodes_array[nodes_index[7]]->FSpring[2] = nodes_array[nodes_index[7]]->FSpring[2]
		- c06*z28;

}

void Element::update_nodes_FPressure() {

	long int node_1, node_2, node_3, node_4;
	double S_x, S_y, S_z;
	double F_x, F_y, F_z;

	//face 1, outward, nodes 3 8 4 7
	node_1 = nodes_index[2];
	node_2 = nodes_index[7];
	node_3 = nodes_index[3];
	node_4 = nodes_index[6];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;

	//face 2, outward, nodes 5 7 6 8
	node_1 = nodes_index[4];
	node_2 = nodes_index[6];
	node_3 = nodes_index[5];
	node_4 = nodes_index[7];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;

	//face 3, outward, nodes 2 7 3 6
	node_1 = nodes_index[1];
	node_2 = nodes_index[6];
	node_3 = nodes_index[2];
	node_4 = nodes_index[5];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;

	//face 4, outward, nodes 1 6 2 5
	node_1 = nodes_index[0];
	node_2 = nodes_index[5];
	node_3 = nodes_index[1];
	node_4 = nodes_index[4];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;

	//face 5, outward, nodes 2 4 1 3
	node_1 = nodes_index[1];
	node_2 = nodes_index[3];
	node_3 = nodes_index[0];
	node_4 = nodes_index[2];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;

	//face 6, outward, nodes 4 5 1 8
	node_1 = nodes_index[3];
	node_2 = nodes_index[4];
	node_3 = nodes_index[0];
	node_4 = nodes_index[7];
	S_x = (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		- (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_y = (-1.0)*(nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[2] - nodes_array[node_4]->position[2])
		+ (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[2] - nodes_array[node_2]->position[2]);
	S_z = (nodes_array[node_1]->position[0] - nodes_array[node_2]->position[0]) * (nodes_array[node_3]->position[1] - nodes_array[node_4]->position[1])
		- (nodes_array[node_3]->position[0] - nodes_array[node_4]->position[0]) * (nodes_array[node_1]->position[1] - nodes_array[node_2]->position[1]);
	F_x = P0*S_x / 8.0;
	F_y = P0*S_y / 8.0;
	F_z = P0*S_z / 8.0;
	nodes_array[node_1]->FPressure[0] = nodes_array[node_1]->FPressure[0] + F_x;
	nodes_array[node_2]->FPressure[0] = nodes_array[node_2]->FPressure[0] + F_x;
	nodes_array[node_3]->FPressure[0] = nodes_array[node_3]->FPressure[0] + F_x;
	nodes_array[node_4]->FPressure[0] = nodes_array[node_4]->FPressure[0] + F_x;
	nodes_array[node_1]->FPressure[1] = nodes_array[node_1]->FPressure[1] + F_y;
	nodes_array[node_2]->FPressure[1] = nodes_array[node_2]->FPressure[1] + F_y;
	nodes_array[node_3]->FPressure[1] = nodes_array[node_3]->FPressure[1] + F_y;
	nodes_array[node_4]->FPressure[1] = nodes_array[node_4]->FPressure[1] + F_y;
	nodes_array[node_1]->FPressure[2] = nodes_array[node_1]->FPressure[2] + F_z;
	nodes_array[node_2]->FPressure[2] = nodes_array[node_2]->FPressure[2] + F_z;
	nodes_array[node_3]->FPressure[2] = nodes_array[node_3]->FPressure[2] + F_z;
	nodes_array[node_4]->FPressure[2] = nodes_array[node_4]->FPressure[2] + F_z;
}

void dump_configuration_vtk(int t_i) {
	long int i, j, k;
	long int node_index;
	long int element_index;

	/////////////////////////////////////////////////////////////////////////////////
	stringstream filename;

	//filename << setw(10) << setfill('0') << simulation_step << ".gel.vtk";
	filename << setw(10) << setfill('0') << t_i << ".gel.vtk";
	ofstream out(filename.str().c_str());
    cout<<"print dump file"<<endl;
	if (!out.is_open()) {
		cout << "Error opening output file " << filename.str().c_str() << endl;
		exit(1);
	}

	out << "# vtk DataFile Version 2.0" << endl;
	out << "Unstructured Grid" << endl;
	out << "ASCII" << endl;
	out << "DATASET UNSTRUCTURED_GRID" << endl;
	out << "POINTS " << N_nodes << " double" << endl;

	for (node_index = 0; node_index < N_nodes; node_index++) {
		out << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[0];
		out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[1];
		out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[2];
		out << endl;
	}

	out << endl;

	out << "CELLS " << N_elements << " " << 9 * N_elements << endl;

	for (element_index = 0; element_index < N_elements; element_index++) {
		out << left << setw(6) << 8;
		for (i = 0; i < 8; i++) {
			out << " " << left << setw(6) << elements[element_index]->nodes_index[i];
		}
		out << endl;
	}

	out << endl;

	out << "CELL_TYPES " << N_elements << endl;

	for (element_index = 0; element_index < N_elements; element_index++) {
		out << left << setw(6) << 12;
		out << endl;
	}

	out << endl;

	out << "CELL_DATA " << N_elements << endl;
	out << "SCALARS cellData double 1" << endl;
	out << "LOOKUP_TABLE default" << endl;
	for (element_index = 0; element_index < N_elements; element_index++) {
		out << right << setw(12) << scientific << setprecision(5) << elements[element_index]->phi;
		out << endl;
	}

	out << endl;

	out.close();

	//////////////////////////////////////////////////////////////////////////////////////
	/*
	stringstream filename2;

	//filename2 << setw(10) << setfill('0') << simulation_step << ".fibers.vtk";
	filename2 << setw(10) << setfill('0') << t_i << ".fibers.vtk";
	ofstream out2(filename2.str().c_str());

	if (!out2.is_open()) {
		cout << "Error opening output file " << filename2.str().c_str() << endl;
		exit(1);
	}

	out2 << "# vtk DataFile Version 2.0" << endl;
	out2 << "polydata" << endl;
	out2 << "ASCII" << endl;
	out2 << "DATASET POLYDATA" << endl;

	long int N_node_on_fibers = 0;

	for (i = 0; i < Nfiber; i++)
		for (j = 0; j < fibers_length[i]; j++) {
			N_node_on_fibers++;
		}

	out2 << "POINTS " << N_node_on_fibers << " double" << endl;

	for (i = 0; i < Nfiber; i++)
		for (j = 0; j < fibers_length[i]; j++) {
			node_index = nodes_on_fibers[i][j];
			out2 << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[0];
			out2 << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[1];
			out2 << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[2];
			out2 << endl;
		}

	out2 << endl;

	long int N_fiber_atom = 0;

	for (i = 0; i < Nfiber; i++) N_fiber_atom += fibers_length[i];

	out2 << "LINES " << Nfiber << " " << N_fiber_atom + Nfiber << endl;

	for (i = 0; i < Nfiber; i++) {
		out2 << left << setw(6) << fibers_length[i];

		for (j = 0; j < fibers_length[i]; j++) {
			out2 << " " << left << setw(6) << nodes[nodes_on_fibers[i][j]]->on_fibers;
		}

		out2 << endl;
	}

	out2 << endl;

	out2.close();
	*/
}

void log_time_phi() {
	int i, j, k;
	long int element_index;

	stringstream filename;

	filename << setw(10) << setfill('0') << simulation_step << ".phi";
	ofstream out(filename.str().c_str());

	if (!out.is_open()) {
		cout << "Error opening output file " << filename.str().c_str() << endl;
		exit(1);
	}

	for (element_index = 0; element_index < N_elements; element_index++) {
		out << left << setw(6) << element_index;
		out << " " << right << setw(12) << scientific << setprecision(5) << elements[element_index]->phi;
		out << endl;
	}

	out.close();
}

void log_node_coordinates(long int node_index) {
	int i, j, k;

	stringstream filename;
	ofstream out;

	filename << node_index <<".coordinates";
	if (simulation_step == 0) {
		out.open(filename.str().c_str());
	}
	else {
		out.open(filename.str().c_str(), ios_base::app);
	}

	if (!out.is_open()) {
		cout << "Error opening output file " << filename.str().c_str() << endl;
		exit(1);
	}

	out << left << setw(6) << simulation_step*h;
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->position[2];
	out << endl;

	out.close();
}

void log_node_velocities(long int node_index) {
	int i, j, k;

	stringstream filename;
	ofstream out;

	filename << node_index << ".velocity";
	if (simulation_step == 0) {
		out.open(filename.str().c_str());
	}
	else {
		out.open(filename.str().c_str(), ios_base::app);
	}

	if (!out.is_open()) {
		cout << "Error opening output file " << filename.str().c_str() << endl;
		exit(1);
	}

	out << left << setw(6) << simulation_step*h;
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->velocity[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->velocity[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->velocity[2];
	out << endl;

	out.close();
}
void log_node_forces(long int node_index) {
	int i, j, k;

	stringstream filename;
	ofstream out;

	filename << node_index << ".forces";
	if (simulation_step == 0) {
		out.open(filename.str().c_str());
	}
	else {
		out.open(filename.str().c_str(), ios_base::app);
	}

	if (!out.is_open()) {
		cout << "Error opening output file " << filename.str().c_str() << endl;
		exit(1);
	}

	out << left << setw(6) << simulation_step*h;
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FSpring[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FSpring[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FSpring[2];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FPressure[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FPressure[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FPressure[2];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCilia[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCilia[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCilia[2];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCiliaBend[0];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCiliaBend[1];
	out << " " << right << setw(12) << scientific << setprecision(5) << nodes[node_index]->FCiliaBend[2];
	out << endl;

	out.close();
}

int load_conf(char *filename) {
	ifstream conf(filename);

	if (!conf.is_open()) {
		cout << "Error opening conf file" << endl;
		exit(1);
	}

	char  buffer[256];
	char *delim = (char *)" ";
	char *tmp;

	int time_written = 0;
	int h_written = 0;
	int dump_written = 0;
	int temperature_written = 0;

	while (!conf.eof())
	{
		conf.getline(buffer, 100);

		if ((strlen(buffer) == 0) || (buffer[0] == '#')) continue;

		tmp = strtok(buffer, delim);

		if (strcmp(tmp, "time") == 0) {
			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout << "conf file error when loading in time: no t_start" <<
					endl;
				exit(1);
			}

			t_start = atof(tmp);
			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout << "conf file error when loading in time: no t_end" <<
					endl;
				exit(1);
			}

			t_end = atof(tmp);
			time_written = 1;
			tmp = strtok(NULL, delim);

			if (tmp != NULL) {
				cout <<
					"conf file error when loading in time: extra argument" <<
					endl;
				exit(1);
			}

			cout << "Time: " << t_start << " ~ " << t_end << endl;
		}
		else if (strcmp(tmp, "h") == 0) {
			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout << "conf file error when loading in h: no argument" <<
					endl;
				exit(1);
			}

			h = atof(tmp);
			h_written = 1;
			tmp = strtok(NULL, delim);

			if (tmp != NULL) {
				cout << "conf file error when loading in h: extra argument" <<
					endl;
				exit(1);
			}
		}
		else if (strcmp(tmp, "dump") == 0) {
			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout <<
					"conf file error when loading in dump: no dump file type" <<
					endl;
				exit(1);
			}

			if (strcmp(tmp, "vtk") == 0) dump_file_type = 1;
			else if (strcmp(tmp, "bin") == 0) dump_file_type = 2;
			else
			{
				cout <<
					"conf file error when loading in dump: dump file type wrong"
					<<
					endl;
				exit(1);
			}

			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout << "conf file error when loading in dump: no dump step" <<
					endl;
				exit(1);
			}

			count_dump_0 = atol(tmp);
			dump_written = 1;
			tmp = strtok(NULL, delim);

			if (tmp != NULL) {
				cout <<
					"conf file error when loading in dump: extra argument" <<
					endl;
				exit(1);
			}
		}
		else if (strcmp(tmp, "temperature") == 0) {
			tmp = strtok(NULL, delim);

			if (tmp == NULL) {
				cout <<
					"conf file error when loading in temperature: no argument"
					<<
					endl;
				exit(1);
			}

			temperature = atof(tmp);
			temperature_written = 1;
			tmp = strtok(NULL, delim);

			if (tmp != NULL) {
				cout <<
					"conf file error when loading in temperature: extra argument"
					<<
					endl;
				exit(1);
			}
		}
		else
		{
			cout << "conf file error: unrecognized line:" << endl;
			cout << buffer << endl;
			//exit(1);
		}
	}

	if (time_written == 0) {
		cout << "conf file error: missing time" << endl;
		exit(1);
	}

	if (h_written == 0) {
		cout << "conf file error: h" << endl;
		exit(1);
	}

	if (dump_written == 0) {
		cout << "conf file error: missing dump" << endl;
		exit(1);
	}

	conf.close();

	return 0;
}
