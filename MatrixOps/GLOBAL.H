/* 
 * Matrix algebra package for linear system solving.
 *
 * Created by: Logan Goldberg
 * Created: 1/8/2017
 * Last Modified: 3/4/2018
 *  
 */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#define EPSILON 1.0e-15
#define BAD_RESULT -1.0e+30
#define MISSING_DATA -1.0e+30

#define SQR(x) ((x) * (x))
#define CUBE(x) (SQR(x) * (x))
#define ONEOVER(x) (1 / (x))
#define RECIPROCAL(x) (ONEOVER(x))

typedef enum { no, yes, dont_care } answer;

#endif
