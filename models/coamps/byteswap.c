/* 
 * This code may (or may not) be part of the COAMPS distribution,
 * So it is not protected by the DART copyright agreement.
 *
 * $Id$
 */

/*******************************
 * MODULE:       coamps_util_mod        
 * AUTHOR:       T. R. Whitcomb
 *               Naval Research Laboratory
 * DART VERSION: Jamaica
 *
 * Routines for in-place byteswapping of FORTRAN arrays
 *******************************/

/******************************
 * BEGIN DEFINITIONS
 ******************************/

#define IN
#define OUT
#define INOUT

#define INLINE

/******************************
 * END DEFINITIONS
 ******************************/

/******************************
 * BEGIN PROTOTYPES
 ******************************/

/* Treat as public routine */
void byteswap_array_( unsigned char*, 
                      int*, 
                      const int*);

/* Treat as private routines */
INLINE void byteswap_value( unsigned char*, 
                            const int);

INLINE void swap_bytes( unsigned char*, 
                        unsigned char*);

/******************************
 * END PROTOTYPES
 ******************************/

/******************************
 * BEGIN ROUTINES
 ******************************/

/**
 * c_byteswap_array
 * ----------------
 * Reverses the byte ordering of an entire array elementwise in 
 * place.  For example:
 *  F90:
 *   real(kind=8), dimension(5) :: test_array
 *   ...
 *   call byteswap_array(test_array, 5, 8)
 *
 * PARAMETERS
 *INOUT array               The array to byteswap
 *   IN array_len           The length of the array (in units of
 *                          element_size)
 *   IN element_size        The size (in bytes) of each element of
 *                          the array
 */
void c_byteswap_array_( INOUT unsigned char *array, 
                        IN    int           *array_len, 
                        IN    const int     *element_size)
{
  /* Need to keep track of where we are both in terms of the actual
   * bytes and in terms of the original elements */
  int ii;
  int element_index;

  for (ii=0; ii < *array_len; ii++)
  {
    element_index = ii * *element_size;
    byteswap_value((array + element_index), *element_size);
  }
}

/**
 * byteswap_value
 * --------------
 * Reverses the byte ordering of a given value in-place
 *
 * PARAMETERS
 *INOUT value               pointer to value to swap bytes in
 *  IN  value_size          The size in bytes of this value 
 *                          (assumed even)
 */
INLINE void byteswap_value( INOUT unsigned char *value, 
                            IN    const int      value_size)
{
  int ii;
  for (ii = 0; ii < (value_size / 2); ii++)
  {
    swap_bytes(value + ii, value + (value_size - ii - 1));
  }
}

/**
 * swap_bytes
 * ----------
 * Switches the value of two parameters
 *
 * PARAMETERS
 *INOUT value1              will hold value2 on exit
 *INOUT value2              will hold value1 on exit      
 */
INLINE void swap_bytes( INOUT unsigned char *value1, 
                        INOUT unsigned char *value2)
{
  unsigned char temp_val;

  temp_val = *value1;
  *value1  = *value2;
  *value2  = temp_val;
}

/******************************
 * END ROUTINES
 ******************************/

/* <next few lines under version control, do not edit>
 * $URL$
 * $Id$
 * $Revision$
 * $Date$
 */
