
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
  cell val;
  int xcoor;
  int ycoor;
  struct linked_cells* next;
} linked_cells;


cell** create_grid(char* seq1, char* seq2){
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    cell** grid = malloc((len1 + 1) * (sizeof(cell*)));
    int i,j;
    for (i = 0; i < (len1 + 1); i++){
        grid[i] = malloc((len2 + 1) * (sizeof(cell)));
    }

    int r = len1 + 1;
    int c = len2 + 1;

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
                grid[i][j].prev = &(grid[i][j]);
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

void print_cell(cell* cell){
  printf("xcoor:%d \t ycoor:%d \t val:%d\n",
          cell->xcoor,cell->ycoor,cell->val);
}

// Write function to free our cells

linked_cells* traceback(cell** grid, int r, int c){
  /* r is the number of rows in the grid
     (not necessarily the number of nucleotides
     in the sequence) and similarly for c */
  printf("Starting traceback\n");
  int x = c;
  //printf("x is initialized\n");
  int y = r;
  //printf("y is initialized\n");
  int step = 0;
  //printf("step is initialized\n");
  
  //printf("A\n");
  linked_cells* rv = malloc(sizeof(linked_cells));
  //printf("B\n");
  linked_cells* curr;
  //printf("C\n");
  curr = rv;
  //printf("D\n");

  while (x != 0 && y != 0){ //terminating condition is we reached the top
    printf("while loop starts\n");
    //first print out (to file?) location
    //print_cell(&(grid[y][x]));
    printf("E\n");
    printf("Row: %d, Col: %d, Val: Step:%d\n", y, x, /*grid[y][x].val ,*/ step);
    
    //fill return value
    curr->xcoor = x;
    curr->ycoor = y;
    curr->val=grid[y][x];
    step++;

    //update parameters
    cell* temp = grid[y][x].prev;
    x = temp->xcoor;
    y = temp->ycoor;
    curr = curr->next = malloc(sizeof(linked_cells));
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
  //f = fopen("/output/seqs.txt", "w+"); 
  // this will need to be changed for multiple files (dynamic assignment)

  linked_cells* temp = path;
  linked_cells* prev = NULL; //this may need to be initialized to something else
  // the idea is we store the previous value to see whether we moved diagonally
  // if we did not move along the proper sequence insert a dash
  
  while(temp != NULL){
    if(temp->xcoor == prev->xcoor)
      fputc('-',f);
    else
      fputc(seq1[temp->xcoor],f);
    prev = temp;
    temp = temp->next;
  } // this prints out seq 1 as aligned

  fputc('\n',f); //this might need to be fprintf - not sure if \n is a char
  temp = path;
  prev = NULL;
  while(temp != NULL){
    if(temp->ycoor == prev->ycoor)
      fputc('-',f);
    else
      fputc(seq2[temp->ycoor],f);
    prev = temp;
    temp = temp->next;
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

void linkedToString(linked_cells* root){
  if (root == NULL)
    printf("end of list\n");
  else {
    printf("(%d,%d) : %d\n",root->xcoor,root->ycoor,root->val.val);
    linkedToString(root->next);
  }
}

int main(int argc, char** argv){
    char* seq1 = "GATTACA";
    char* seq2 = "GCATGCU";
    cell** grid = create_grid(seq1, seq2);
    print_grid(grid, 8, 8);
    printf("does this work?\n");
    linked_cells* l1 = traceback(grid, 8, 8);
    printf("Passed traceback\n");
    linked_cells* l1_r = reverse(l1);
    printf("passed reverse\n");
    FILE* f = fopen("/output/seqs.txt", "w+");
    printSeqToFile(f,l1_r,seq1,seq2);
    fclose(f);
    free_grid(grid, 8, 8);
    free_linked_cells(l1);
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
