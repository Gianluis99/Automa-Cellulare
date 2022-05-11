#include <iostream>
#include<mpi.h>
#include<cstdio>
#include<cmath>
#include<vector>
#include<fstream>
#include <allegro5/allegro.h>
#include<allegro5/allegro_primitives.h>
#include<allegro5/allegro_native_dialog.h>
#include<allegro5/allegro_font.h>
#include<allegro5/allegro_ttf.h>	

using namespace std;

//COSTANTI
const  int EMPTY=0;
const  int GAS=1; //blue
const  int ICE=2; //yellow
const  int FIRE=3; //red
const  int BURNT=4; 
const  int N=100;
const int passi=25;
int arrayRead[(N*N)];



/*tre strati:vuoto, gas, e congelato. il gas avanza automaticamente in celle libere, se tocca una cella ghiacciata allora si ghiaccia.*/


struct Coor{
int i;
int j;

Coor(int i_,int j_){
i=i_;
j=j_;
}

};
void check(int i,int j,int *arrayReadl,int *arrayWritel){
 
 
  if(arrayReadl[i*N+j]==GAS){
  	

  //fire
  
  	bool fire=false;
   if ( arrayReadl[(i*N)+(j-1)]==3 )
 	{fire=true;}
 	
if ( arrayReadl[(i*N)+(j+1)]==3 )
	 {fire=true;}
	 
if ( arrayReadl[((i+1)*N)+j]==3 )
	 {fire=true;}
	 
if ( arrayReadl[((i+1)*N)+(j+1)]==3)
	{fire=true;}
	
if( arrayReadl[((i+1)*N)+(j-1)]==3)
	{fire=true;} 
	
if ( arrayReadl[((i-1)*N)+j]==3 )
	{fire=true;}
	
if ( arrayReadl[((i-1)*N)+(j+1)]==3)
	{fire=true;}
	
if ( arrayReadl[((i-1)*N)+(j-1)]==3)
	{fire=true;}
  
  
  if(fire)
  {
	arrayWritel[((i-1)*N)+(j-1)]=4;
	arrayWritel[((i-1)*N)+(j+1)]=4;
	arrayWritel[((i-1)*N)+j]=4;
	arrayWritel[((i+1)*N)+(j-1)]=4;
	arrayWritel[((i+1)*N)+(j+1)]=4;
	arrayWritel[((i+1)*N)+j]=4;
	arrayWritel[(i*N)+(j+1)]=4;
	arrayWritel[(i*N)+(j-1)]=4;
	arrayWritel[(i*N)+(j-1)]=4;
	arrayWritel[i*N+j]=4;
	return;

  }
  
     if ( arrayReadl[(i*N)+(j-1)]==2 )
 	{arrayWritel[i*N+j]=2;  return;}
 	
if ( arrayReadl[(i*N)+(j+1)]==2 )
	 {arrayWritel[i*N+j]=2;  return;}
	 
if ( arrayReadl[((i+1)*N)+j]==2 )
	 {arrayWritel[i*N+j]=2;  return;}
	 
if ( arrayReadl[((i+1)*N)+(j+1)]==2)
	{arrayWritel[i*N+j]=2;	return;}
	
if( arrayReadl[((i+1)*N)+(j-1)]==2)
	{arrayWritel[i*N+j]=2;  return;} 
	
if ( arrayReadl[((i-1)*N)+j]==2 )
	{arrayWritel[i*N+j]=2;  return;}
	
if ( arrayReadl[((i-1)*N)+(j+1)]==2)
	{arrayWritel[i*N+j]=2;  return;}
	
if ( arrayReadl[((i-1)*N)+(j-1)]==2)
	{arrayWritel[i*N+j]=2;  return;}
  
  
    	  
    	  vector<Coor>v;

  	 if ( arrayReadl[(i*N)+(j-1)]==0 )
 	{   Coor c(i,j-1);
 	    v.push_back(c);
 		}

	if ( arrayReadl[(i*N)+(j+1)]==0 )
	{	Coor c(i,j+1);
 		v.push_back(c);
 	}

	 

	 if ( arrayReadl[((i+1)*N)+j]==0)
	{	
	Coor c(i+1,j);
 		v.push_back(c);
 	}

	 

	 if ( arrayReadl[((i-1)*N)+j]==0 )
	{	
		Coor c(i-1,j);
 		v.push_back(c);
 	}

	
	if(!v.empty()){
	 arrayWritel[i*N+j]=0; 
	 
	 int num=((i*j)+rand())%v.size();
    	   arrayWritel[(v[num].i)*N+(v[num].j)]=1; 
    	  
    	 }
	else
	arrayWritel[i*N+j]=GAS;

	return;
 }

 else if(arrayReadl[i*N+j]==ICE)arrayWritel[i*N+j]=ICE;
 else if(arrayReadl[i*N+j]==FIRE){
 if ( arrayReadl[(i*N)+(j-1)]==4 )
 	{arrayWritel[i*N+j]=4;  return;}
 	
if ( arrayReadl[(i*N)+(j+1)]==4 )
	 {arrayWritel[i*N+j]=4;  return;}
	 
if ( arrayReadl[((i+1)*N)+j]==4 )
	 {arrayWritel[i*N+j]=4;  return;}
	 
if ( arrayReadl[((i+1)*N)+(j+1)]==4)
	{arrayWritel[i*N+j]=4;	return;}
	
if( arrayReadl[((i+1)*N)+(j-1)]==4)
	{arrayWritel[i*N+j]=4;  return;} 
	
if ( arrayReadl[((i-1)*N)+j]==4 )
	{arrayWritel[i*N+j]=4;  return;}
	
if ( arrayReadl[((i-1)*N)+(j+1)]==4)
	{arrayWritel[i*N+j]=4;  return;}
	
if ( arrayReadl[((i-1)*N)+(j-1)]==4)
	{arrayWritel[i*N+j]=4;  return;}
	
	arrayWritel[i*N+j]=FIRE;
 }
 else if(arrayReadl[i*N+j]==BURNT)arrayWritel[i*N+j]=EMPTY;
 	
 	
}



void copy( int *arrayRead,int *arrayWrite,int x){
	for(int i=N;i<x+N;i++)
		arrayRead[i]=arrayWrite[i];


}



void inizializza(){

	srand (time(NULL));
    	int num;

	     
	  for(int i=0;i<N;i++){
	  	for(int j=0;j<N;j++){
	  	  if((i==0 || j==0) || i==N-1||j==(N-1) || (i==N/2 &&j==N/2) || (i==N/2+1 &&j==N/2+1)|| (i==N/2+2 &&j==N/2+1)
	  	   ){
	  	         arrayRead[i*N+j] = ICE;
	  	  }
	  	  else{	
	  	  		
	    num=rand()%30;
           
            if(num >=2 && num<=6){
                arrayRead[i*N+j] = GAS;
                }
                else if(num ==0  )
                  arrayRead[i*N+j] =FIRE;
                 else if(num >=7 && num<8)
                    arrayRead[i*N+j] =ICE;
                  else
                   arrayRead[i*N+j] = EMPTY;
	    }
	   }
	 }

}

void printChunkRank(int * localArray,int x){

for(int i=N;i<x+N;i++)
    cout<<localArray[i];

}
   

void draw(){

 for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            if (arrayRead[x*N+y] == 1) {
            	al_draw_filled_rectangle(x * 4 + 4, y * 4+4, x * 4, y * 4 ,al_map_rgb(0,0,255));

            }
            else if (arrayRead[x*N+y] == 2) {
                 al_draw_filled_rectangle(x * 4 + 4, y * 4+4, x *4, y * 4 ,al_map_rgb(255, 211, 25));
            }
             else if (arrayRead[x*N+y] == 0) {
                 al_draw_filled_rectangle(x * 4 + 4, y * 4+4, x * 4, y * 4 ,al_map_rgb(0, 0, 0));
            }
            else if (arrayRead[x*N+y] == 3) {
                 al_draw_filled_rectangle(x * 4 + 4, y * 4+4, x * 4, y * 4 ,al_map_rgb(255, 0, 0));
            }
             else if (arrayRead[x*N+y] == 4) {
                 al_draw_filled_rectangle(x * 4 + 4, y * 4+4, x * 4, y * 4 ,al_map_rgb(108, 14, 89));
            }
        }
    }
    	
al_flip_display();
al_rest(0.6);
}





//MPI

int main(int argc, char *argv[]){

int *localArrayRead;
int *localArrayWrite;
int numRows,numCol;

//topologia virtuale
 MPI_Comm top; // nuovo comunicatore
  int dim[2],period[2],reorder; //dati comunicatore
   int coord[2],id;

 //tipo derivato
 MPI_Datatype row; 
 
int e,nump;
int x;
int myid;
int dest,source;

int up,down;
double startTime,endTime,time;

MPI_Request requests[4];
MPI_Status statuses[4];

ALLEGRO_DISPLAY *display;

 
e=MPI_Init( &argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD,&nump);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);

//***topologia******
dim[0]=1; dim[1]=nump; //dimensione 0=sono colonne, dimensione 1=sono righe

period[0]=false; period[1]=false;
reorder=true;
MPI_Cart_create(MPI_COMM_WORLD,2,dim,period,reorder,&top);


//inizializzazione matrice e display da parte di un solo thread
if(myid==0){
	
if(!al_init())
{
	cerr<<"erorr";
	return -1;
}

al_set_new_display_flags(ALLEGRO_NOFRAME);
display=al_create_display(400,400);
al_set_window_position(display,300,200);
al_set_window_title(display,"Automa SOOT");
 al_init_primitives_addon();
al_clear_to_color(al_map_rgb(0,0,0));

inizializza();
draw();
}	
x=(N*N)/nump;  //andiamo a ottenere le porzioni da dare ai thread
localArrayRead=new int [x+(N*2)]; //aggiungiamo orlatura  ovvero due righe da N	
localArrayWrite=new int [x+(N*2)]; //aggiungiamo orlatura 	
 numRows=N/nump;

 
//dopo l'inizializzazione  e allocazione mettiamo una barriera
MPI_Barrier(top); 
//startTime=MPI_Wtime();

//creiamo nuovo tipo derivato ovvero la riga
MPI_Type_contiguous(N,MPI_INT,&row);
MPI_Type_commit(&row);



MPI_Scatter(&arrayRead[0],x,MPI_INT, &localArrayRead[N],x,MPI_INT,0,top); //Ricevono a partire dalla seconda riga
    	
  //Cerchiamo i vicini a cui spedire
    MPI_Cart_shift(top,1,1,&up,&down); 
  // printf("P:%d i vicini sono:  up:%d  down:%d \n",myid,up,down );
	
	for(int p=0;p<passi;p++){
	    //se il thread ha l'up gli invia la riga
  	 	 MPI_Isend(&localArrayRead[N],1,row,up,10,top,&requests[0]); //gli invio la mia 1 riga
  	 	 MPI_Irecv(&localArrayRead[0],N,MPI_INT,up,MPI_ANY_TAG,top,&requests[1]);
  	 	 
            //se il thread ha il down
  	 	 MPI_Isend(&localArrayRead[(numRows)*N],1,row,down,11,top,&requests[2]); //spedisce la penultima la 25
         	 MPI_Irecv(&localArrayRead[(numRows+1)*N],N,MPI_INT,down,MPI_ANY_TAG,top,&requests[3]); //nell'ultima riga ovvero 26
    		
    		 
		//inizialmente non consideriamo la prima riga, non consideriamo l'ultima
		for(int i=2;i<numRows;i++){ 
		   for(int j=0;j<N;j++){
			check(i,j,localArrayRead,localArrayWrite);
			}
		}
		
		
         MPI_Waitall(4,requests,statuses);
	         
	//calcola prima riga
  		   for(int i=1;i<2;i++){ 
		   for(int j=0;j<N;j++){
			check(i,j,localArrayRead,localArrayWrite);
			}
		  }
	//calcola l'ultima riga
            for(int i=numRows;i<numRows+1;i++){ //si calcola la 25esima 
		   for(int j=0;j<N;j++){
			check(i,j,localArrayRead,localArrayWrite);
			}
		  }
		  
  		copy(localArrayRead,localArrayWrite,x);
		
		
		MPI_Gather(&localArrayWrite[N],x,MPI_INT,&arrayRead,x,MPI_INT,0,top);
		
		if(myid==0){
		draw();
		}
			
	
}
/*
endTime=MPI_Wtime();
time=endTime-startTime;
 
  */
if(myid==0){
	al_destroy_display(display);
}

delete[] localArrayRead;
delete[] localArrayWrite;
MPI_Type_free(&row);//andiamo a deallocare il tipo derivato
e=MPI_Finalize(); //end


}







