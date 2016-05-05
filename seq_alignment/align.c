
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
    for (int i = 0; i < (len1 + 1); i++){
        grid[i] = malloc((len2 + 1) * (sizeof(cell)));
    }

    int r = len1 + 1;
    int c = len2 + 1;

    for (int i = 0; i < r; i++){        
        grid[i][0].val = GAP * i;
        if (i != 0){
            grid[i][0].prev = &(grid[i-1][0]);
        } 
        for (int j = 1; j < c; j++){
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
    for (int i = 0;i < r; i++){
        for (int j = 0; j < c; j++){
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

  int x = c;
  int y = r;
  int step = 0;
  //adopting convention that the x coordinate is column index
  
  linked_cells* rv = malloc(sizeof(linked_cells));
  linked_cells* curr;
  curr = rv;
  while (x != 0 && y != 0){ //terminating condition is we reached the top

    //first print out (to file?) location
    printf("Row: %d, Col: %d, Val:%d, Step:%d\n", y, x, grid[y][x].val, step);
    
    //fill return value
    curr->xcoor = x;
    curr->ycoor = y;
    curr->val=grid[y][x];
    step++;

    //update parameters
    cell* temp = grid[y][x].prev;
    x = temp->xcoor;
    y = temp->ycoor;
    curr = curr->next;
  }
  return rv;
}

void free_linked_cell(linked_cell* tbd){
  if (tbd->next != NULL)
    free_linked_cell(tbd->next);
  free(tbd);
}

void free_grid(cell** grid,int r, int c){
  for(int i = 0; i < r; i++){
  /*  for(int j = 0; j < c; j++)
      free(grid[i][j]); //not sure if we need this
  */
    free(grid[i]);
  }
  free(grid)
}

//Need to update cell struct and be wary of convention with xcoor/ycoor


int main(int argc, char** argv){
    char* seq1 = "GATTACA";
    char* seq2 = "GCATGCU";
    cell** grid = create_grid(seq1, seq2);
    print_grid(grid, 8, 8);
    free_grid(grid);

    return 0;
} 

// Need to write len_linked_cells
// Need to write linked_cells_invert
// Might want print_linked_cells
// Need to write to file from traceback
// Need to look at passing around lengths
// Decide on input and output format



// Right now we're checking the length of the sequences

// Write the traceback function

// Output the aligned sequences

// Userfriendly

// More nuanced scoring 