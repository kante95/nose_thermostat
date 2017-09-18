#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>


#define N 125

double myrandom(){
	return ((double)rand()/(double)(RAND_MAX));
}

double normrand(double std){
	return sqrt((-2*(std*std)*log(1-myrandom()))*cos(3.1415*2*myrandom()));
} 

void printmatrix(double matrix[N][3]){
	int i;
	for(i=0;i<N;i++){
		printf("Particella %d: [%lf] [%lf] [%lf]\n",i+1,matrix[i][0],matrix[i][1],matrix[i][2]);

	}

}

void copy_matrix(double source[N][3],double destination[N][3]){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<3;j++){
			destination[i][j] = source[i][j];
		}
	}
}


void initialize(double L,double positions[N][3],double velocities[N][3]){
    float i = 5;
    int index = 0;
    float j,k,m;
    for(j =0;j<i;j++){
    	for(k =0;k<i;k++){
    		for(m =0;m<i;m++){
    			positions[index][0]= j*(L/i) - L/2.0;
                positions[index][1] = k*(L/i) - L/2.0;
                positions[index][2] = m*(L/i) - L/2.0;
                velocities[index][0]= normrand(1);
                velocities[index][1] = normrand(1);
                velocities[index][2] = normrand(1);
                //printf("x: %lf y: %lf z: %lf",positions[index][0],positions[index][1],positions[index][2]);
                index++;
			}
    	}
    }
}

double * lennard_jones_derivative(double x,double y,double z){
    double den =pow(x,2)+pow(y,2) +pow(z,2);
    double static results[3];
    results[0] = -4*(-12*x*pow(den,-7) + 6*x*pow(den,-4));
    results[1] = -4*(-12*y*pow(den,-7) + 6*y*pow(den,-4));
    results[2] = -4*(-12*z*pow(den,-7) + 6*z*pow(den,-4));
    return results;
}


double lennard_jones(double r){
	double r2 = r*r;
	double r6 = r2*r2*r2;
	double r12 = r6*r6;
	return 4*(1.0/r12-1.0/r6);
}

double calculate_potential(double positions[N][3],double L,double acceleration[N][3]){

    double potential = 0;
    double r;
    int i,j,x;
    double d,a,b,c;
    double *force;
    for(i=0;i<N;i++){
    	for(j=i+1;j<N;j++){
    		r=0;
    		for(x = 0; x<3;x++){
    			d = positions[i][x] - positions[j][x]-L*rint((positions[i][x] - positions[j][x])/L);
    			r+=d*d;
    		}
    		r = sqrt(r);
    		if(r<L/2.0 && r!=0.0){
    			potential+=2*lennard_jones(r);
    			a = positions[i][0] - positions[j][0]-L*rint((positions[i][0] - positions[j][0])/L);
    			b = positions[i][1] - positions[j][1]-L*rint((positions[i][1] - positions[j][1])/L);
    			c = positions[i][2] - positions[j][2]-L*rint((positions[i][2] - positions[j][2])/L);
    			force =  lennard_jones_derivative(a,b,c);
    			for(x = 0; x<3;x++){
    				acceleration[i][x] = force[x];
    				acceleration[j][x] = -force[x];
    			}
    		}
    	}
    }

    return potential;
}

int simulation(double density,double temperature,double time, double step){
    double V = N/density;
    double L = pow(V,1.0/3.0);

    double positions[N][3];
    double velocities[N][3];

    double current_accelerations[N][3];
    double new_accelerations[N][3];

    double potential,energy;

    int t,i,j;
   
    printf("Inizializzo lo stato iniziale....\n");
    initialize(L,positions,velocities);
    potential = calculate_potential(positions,L,current_accelerations);
	energy = potential;
	for(i=0;i<N;i++){
			for(j=0;j<3;j++){
				energy += 0.5* pow(velocities[i][j],2);
		}
	}

    
    for(t = 0; t<time;time+=step){
        //print("$positions\n")
        //print("Simulando secondo "+str(t)+" s")
        //calcolo le nuove posizioni
        //print("$shape(positions) $shape(velocities)")
        for(i=0;i<N;i++){
			for(j=0;j<3;j++){
				positions[i][j] = positions[i][j] + velocities[i][j]*step + 0.5*current_accelerations[i][j]*pow(step,2);
        		positions[i][j] = positions[i][j] - L*rint(positions[i][j]/L);
			}
		}
        
        //accelerazioni dovute alle nuove posizioni
        potential = calculate_potential(positions,L,new_accelerations);
        
        for(i=0;i<N;i++){
			for(j=0;j<3;j++){
        		velocities[i][j] = velocities[i][j] + 0.5*(current_accelerations[i][j] + new_accelerations[i][j])*step;
    		}
    	}
        //aggiorno le accelerazioni per il prossimo ciclo
        copy_matrix(new_accelerations,current_accelerations);
 
        energy+= potential;
        for(i=0;i<N;i++){
			for(j=0;j<3;j++){
				energy += 0.5* pow(velocities[i][j],2);
			}
		}
        printf("Energy at second %f: %f\n",t,energy);
    }
    return 0;
}

int main(){
	simulation(0.01,1,0.1,0.00001);
    
    return 0;
}