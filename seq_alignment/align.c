
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int MATCH = 1;
int MISMATCH = -1;
int GAP = -1;

typedef struct cell{
    int val;
    int xcoor;
    int ycoor;
    struct cell* prev;
} cell;

typedef struct linked_cells{
  int xcoor;
  int ycoor;
  struct linked_cells* next;
} linked_cells;

void print_cell(cell* cell){
  printf("xcoor:%d \t ycoor:%d \t val:%d\n",
          cell->xcoor,cell->ycoor,cell->val);
}

void linkedToString(linked_cells* root){
  if (root == NULL)
    printf("end of list\n");
  else {
    printf("(%d,%d)\n",root->xcoor,root->ycoor);
    linkedToString(root->next);
  }
}
/*
cell** dyn_grid_gen(char* seq1, char* seq2){
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    cell** grid = (cell**)malloc((len1 + 1) * (sizeof(cell*)));
    int i,j;
    for (i = 0; i < (len1 + 1); i++){
        grid[i] = (cell*)malloc((len2 + 1) * (sizeof(cell)));
    }

    int r = len1 + 1;
    int c = len2 + 1;

}
*/


cell** create_grid(char* seq1, char* seq2, bool dynamic){
    // lets say we add a boolean 'dynamic' that determines the gap/mismatch behavior
    // how do we implement efficiently?
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    cell** grid = (cell**)malloc((len1 + 1) * (sizeof(cell*)));
    int i,j;
    for (i = 0; i < (len1 + 1); i++){
        grid[i] = (cell*)malloc((len2 + 1) * (sizeof(cell)));
    }

    int r = len1 + 1;
    int c = len2 + 1;

    if (dynamic){

    } else { 

        for (i = 0; i < r; i++){    
            grid[i][0].val = GAP * i;
            if (i != 0){
                grid[i][0].prev = &(grid[i-1][0]);
            } 
            for (j = 1; j < c; j++){
                if (i == 0){
                    grid[i][j].val = GAP * j;
                    if (j != 0){
                        grid[i][j].prev = &(grid[i][j-1]);
                        }
                } else {
                    char nucl1 = seq1[i-1];
                    char nucl2 = seq2[j-1];
                    int diagonal;
                    if (nucl1 == nucl2){
                        diagonal = MATCH;
                    } else {
                        diagonal = MISMATCH;
                    }
                    int max = diagonal + grid[i-1][j-1].val;
                    grid[i][j].prev = &(grid[i-1][j-1]);
                    if (GAP + grid[i-1][j].val > max){
                        max = GAP + grid[i-1][j].val;
                        grid[i][j].prev = &(grid[i-1][j]);
                    }
                    if (GAP + grid[i][j-1].val > max){
                        max = GAP + grid[i][j-1].val;
                        grid[i][j].prev = &(grid[i][j-1]);
                    }
                    grid[i][j].val = max;
                    grid[i][j].xcoor = j;
                    grid[i][j].ycoor = i;
                }
            }
        }
    }  
    return grid;
}

void print_grid(cell** grid, int r, int c){
    int i,j;
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            printf("%d\t", grid[i][j].val);
        }
        printf("\n");
    }
}


// Write function to free our cells

linked_cells* traceback(cell** grid, int r, int c){
  /* r is the number of rows in the grid
     (not necessarily the number of nucleotides
     in the sequence) and similarly for c */
  int x = c - 1;
  int y = r - 1;
  int step = 0;
  
  linked_cells* rv = (linked_cells*)malloc(sizeof(linked_cells));
  //printf("B\n");
  linked_cells* curr;
  //printf("C\n");
  curr = rv;

  while (x != 0 && y != 0){ //terminating condition is we reached the top
    
    //fill return value
    curr->xcoor = x;
    curr->ycoor = y;
    step++;

    //update parameters
    cell* temp = grid[y][x].prev;
    x = temp->xcoor;
    y = temp->ycoor;
    curr = curr->next = (linked_cells*)malloc(sizeof(linked_cells));
  }
  return rv;
}

linked_cells* reverse(linked_cells* root){
  linked_cells* temp;
  linked_cells* previous = NULL;
  while(root != NULL) {
    temp = root->next;
    root->next = previous;
    previous = root;
    root = temp;
  }
  return previous;
}

void printSeqToFile(FILE* f, linked_cells* path, char* seq1, char* seq2){
    linked_cells* root = path;
    linked_cells* temp = root->next;

    while(temp != NULL){
        if(root->xcoor == temp->xcoor){
            fputc('-', f);
        } else {
            fputc(seq2[(temp->xcoor) - 1], f);
        }
        root = temp;
        temp = root->next;
    } // this prints out seq 1 as aligned

    fputc('\n',f); //this might need to be fprintf - not sure if \n is a char
    root = path;
    temp = root->next;
    while(temp != NULL){
        if(root->ycoor == temp->ycoor){
            fputc('-',f);
        } else {
            fputc(seq1[(temp->ycoor) - 1],f);
        }
        root = temp;
        temp = root->next;
    } // this prints out seq 2 as aligned
    return;
}


void free_linked_cells(linked_cells* tbd){
  if (tbd->next != NULL)
    free_linked_cells(tbd->next);
  free(tbd);
}

void free_grid(cell** grid, int r, int c){
  int i;
  for(i = 0; i < r; i++){
  /*  for(int j = 0; j < c; j++)
      free(grid[i][j]); //not sure if we need this
  */
   free(grid[i]);
  }
  free(grid);
}

//Need to update cell struct and be wary of convention with xcoor/ycoor

int len_linked_cells(linked_cells* root){
  if (root == NULL)
    return 0;
  else
    return 1 + len_linked_cells(root->next);
}


int main(int argc, char** argv){
    char* seq1 = "GATTACA";
    char* seq2 = "GCATGCU";
    cell** grid = create_grid(seq1, seq2);
    print_grid(grid, 8, 8);
    linked_cells* l1 = traceback(grid, 8, 8);
    linked_cells* l1_r = reverse(l1);
    linkedToString(l1_r);
    FILE* f = fopen("output/seqs.txt", "w+");
    printSeqToFile(f,l1_r,seq1,seq2);
    fclose(f);
    free_grid(grid, 8, 8);
    //free_linked_cells(l1); this is a redundant free statement
    free_linked_cells(l1_r);
    return 0;
} 

/* Need to write len_linked_cells          --done
   Need to write linked_cells_invert       --done
   Might want print_linked_cells           --done
   Need to write to file from traceback    --done
   Need to look at passing around lengths  --
   Decide on input and output format       --


   Right now we're checking the length of the sequences
      - This might need to be made into a global variable
      - a lot of functions need it and future ones might too

   *DONE* Write the traceback function

   *DONE NEEDS TESTING* Output the aligned sequences

   Userfriendly Question, who is this code for?

   More nuanced scoring 
*/
