//
//  prettyFile.c
//  
//
//  Created by Stephen Dellinger on 4/8/13.
//
//

#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[]){
    FILE *in;
    
    //ask for file name
    printf("\nEnter the name of the input file > ");
    ni = scanf("%s",infile);
    
    //open the in file for use
    in = fopen(infile,"r");
    
    //error checking for good filenames
    if((in = fopen(infile,"r")) == NULL) {
        printf("\nCannot open file for input\n");
        exit(1);
    }
    
    while (fscanf(in,) != EOF) {}
}