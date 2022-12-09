
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// interface graphique
#include "graph.h"

#define CARRE 3
#define INTER 0
#define HM 240
#define LM 400
unsigned long noir,vert;

// nombre max d'itérations
#define ITER 10000

/* longueur cycle recherche - 1 */
#define LONGCYCLE 11

#define DIFFTEMPS(a,b) (((b).tv_sec - (a).tv_sec) + ((b).tv_usec - (a).tv_usec)/1000000.)

typedef char tab[HM][LM];


void init( tab );
void calcnouv( tab, tab );
void affiche( tab );

//variable globale pour éviter le stack overflow
tab tt[LONGCYCLE];

int main()
{
	int i;
	struct timeval tv_init, tv_end;

	/* init interface graphique */
	initgraph(LM*(CARRE+INTER), HM*(CARRE+INTER));
	vert = newcolor(0.,0.6,0.);
	noir = newcolor(0.,0.,0.);

	/* initialisation du premier tableau */
	init( tt[0] );

	/* mesure temps */
	gettimeofday( &tv_init, NULL);

	/* itérations */
	for( i=0 ; i<ITER ; i++ )
	{
		affiche( tt[i%LONGCYCLE] );
		calcnouv( tt[i%LONGCYCLE], tt[(i+1)%LONGCYCLE] );

		// à commenter pour ne pas attendre entre deux itérations
		//if( waitgraph() == XK_q )
			// l'utilisateur a appuyé la touche 'q'
			//break;

	}

	gettimeofday( &tv_end, NULL);
	/* temps d'exécution */
	printf( "Temps : %lfs.\n", DIFFTEMPS(tv_init,tv_end) );

	/* fermeture interface graphique */
	closegraph();
	return( 0 );
}

// initialisation d'un tableau qui évolue de manière intéressante
void init( tab t )
{
	int i,j;

/*	srand(time(0)); */
	for( i=0 ; i<HM ; i++ )
		for( j=0 ; j<LM ; j++ )
		{
//			t[i][j] = rand()%2;
//			t[i][j] = ((i+j)%3==0)?1:0;
//			t[i][j] = (i==0||j==0||i==HM-1||j==LM-1)?0:1;
			t[i][j] = 0;
		}

/*	// un marcheur en haut à gauche, qui va percuter ...
	t[10][10] = 1;
	t[10][11] = 1;
	t[10][12] = 1;
	t[9][12] = 1;
	t[8][11] = 1;

	// ... une structure qui se stabilise rapidement vers (50,50)
	t[55][50] = 1;
	t[54][51] = 1;
	t[54][52] = 1;
	t[55][53] = 1;
	t[56][50] = 1;
	t[56][51] = 1;
	t[56][52] = 1;
*/

	// un "mathusalhem" : le gland.
	t[150][100] = 1;
	t[151][100] = 1;
	t[151][102] = 1;
	t[153][101] = 1;
	t[154][100] = 1;
	t[155][100] = 1;
	t[156][100] = 1;

	// en bas : une structure avec un cycle de 3 (une croix)
	int croix[][2]= { {0,3}, {1,3}, {1,2}, {1,1} };	// un 8ième de croix
	for( i=0 ; i<4 ; i++ ) // les 4 coordonnées du tableau
		for( j=0 ; j<4; j++ ) // multipliées 4*2 fois
		{
			t[230+((j%2)?(1+croix[i][0]):(-croix[i][0]))][200+((j<2)?(1+croix[i][1]):(-croix[i][1]))] = 1;
			t[230+((j%2)?(1+croix[i][1]):(-croix[i][1]))][200+((j<2)?(1+croix[i][0]):(-croix[i][0]))] = 1;
		}
	// à droite : un octogone à cycle 5
	int octo[][2]= { {0,3}, {1,2} };	// un 8ième de croix
	for( i=0 ; i<2 ; i++ ) // les 4 coordonnées du tableau
		for( j=0 ; j<4; j++ ) // multipliées 4*2 fois
		{
			t[150+((j%2)?(1+octo[i][0]):(-octo[i][0]))][300+((j<2)?(1+octo[i][1]):(-octo[i][1]))] = 1;
			t[150+((j%2)?(1+octo[i][1]):(-octo[i][1]))][300+((j<2)?(1+octo[i][0]):(-octo[i][0]))] = 1;
		}

}

// calcule le nombre de voisins de la case (i,j) du tableau t.
int inline nbvois( tab t, int i, int j )
{
	int n=0;
	if( i>0 )
	{	/* i-1 */
		if( j>0 )
			if( t[i-1][j-1] )
				n++;
		if( t[i-1][j] )
				n++;
		if( j<LM-1 )
			if( t[i-1][j+1] )
				n++;
	}

		/* i */
	if( j>0 )
		if( t[i][j-1] )
			n++;
	if( j<LM-1 )
		if( t[i][j+1] )
			n++;

	if( i<HM-1 )
	{	/* i+1 */
		if( j>0 )
			if( t[i+1][j-1] )
				n++;
		if( t[i+1][j] )
				n++;
		if( j<LM-1 )
			if( t[i+1][j+1] )
				n++;
	}
	return( n );
}

// calcule un nouveau tableau (n) en fonction d'un ancien (t)
void calcnouv( tab t, tab n )
{
	int i,j, v;
	for(  i=0 ; i<HM ; i++ )
	{
		for( j=0 ; j<LM ; j++ )
		{
			v = nbvois( t, i, j );
			if( v==3 )
				n[i][j] = (t[i][j])?1:2;
			else if( v==2 )
				n[i][j] = (t[i][j])?1:0;
			else
				n[i][j] = 0;
		}
	}
}

// affiche un tableau t
void affiche( tab t )
{
	int i,j,k;
	cleargraph();
	for( i=0 ; i<HM ; i++ )
	{
		for( j=0 ; j<LM ; j++ )
		{
			if( t[i][j]==1 )
				setcolor( noir );
			else if( t[i][j]==2 )
				setcolor( vert );
			if( t[i][j] )
				for( k=0 ; k<CARRE ; ++k )
					line( j*(CARRE+INTER)+k,i*(CARRE+INTER),
					      j*(CARRE+INTER)+k,i*(CARRE+INTER)+CARRE-1 );
		}
	}
	refresh();
}
