/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 *
 * DART $Id$
 */

/* convert gregorian days/seconds to and from year/month/day/hr/min/sec
 *
 * usage:  gregorian_time  days  seconds
 *    or:  gregorian_time  year month day hour minute seconds
 *
 * nsc 12sep2007
 *     17oct2007
 *     20nov2007
 */

#include <stdio.h>
#include <stdlib.h>

/* the start of gregorian numbered days */
#define BASE_YEAR 1601

enum todo { TO_YMD, TO_DS };
enum err  { ERROR, OK };

/* 13th month is leap feb, and use month number directly so
 * add a dummy pad for month 0.
 */
int days_per_month[14] = 
{ 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 29 };

void convert(int direction, int *days, int *secs, 
                            int *year, int *month, int *day,
                            int *hour, int *min, int *sec);
enum err vet_time(int days, int secs);
enum err vet_date(int year, int month, int day, int hour, int min, int sec);

main(int argc, char **argv)
{
    int days, secs, year, month, day, hour, min, sec;
    enum err error;
    enum todo direction;

    /* going from days/sec to y/m/d/h/m/s */
    if (argc == 3) {
        direction = TO_YMD;
        days = atoi(argv[1]);
        secs = atoi(argv[2]);
    } else if (argc == 7) {
        direction = TO_DS;
        year  = atoi(argv[1]);
        month = atoi(argv[2]);
        day   = atoi(argv[3]);
        hour  = atoi(argv[4]);
        min   = atoi(argv[5]);
        sec   = atoi(argv[6]);
    } else {
        fprintf(stderr, "Converts time between gregorian days/seconds (days since 1601)\n");
        fprintf(stderr, "  and YYYY/MM/DD HH:MM:SS (and back)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s days seconds\n", argv[0]);
        fprintf(stderr, "   or  %s year month day hour min sec\n", argv[0]);
        fprintf(stderr, "\n");
        exit (-1);
    }


    if (direction == TO_YMD)
        error = vet_time(days, secs);
    else
        error = vet_date(year, month, day, hour, min, sec);

    if (error != OK)
        exit(-2);


    convert(direction, &days, &secs, &year, &month, &day, 
                       &hour, &min, &sec);

    printf("%6d %5d == %4d/%02d/%02d %02d:%02d:%02d\n", days, secs,
            year, month, day, hour, min, sec);

    exit (0);
} 

/* make sure days/secs have valid range */
enum err vet_time(int days, int secs) {
    /* assume days can be anything, relative to 1601, but must be positive
     * because the convert routine will not go backwards.
     */
    if (days < 0) {
        fprintf(stderr, "error: days must be non-negative\n");
        return ERROR;
    } 

    /* seconds per day */
    if (secs < 0 || secs >= 86400) {
        fprintf(stderr, "error: seconds must be between 0 and 86399\n");
        return ERROR;
    } 
    return OK;
}

/* make sure y/m/d h:m:s is valid */
enum err vet_date(int year, int month, int day, int hour, int min, int sec) {

    int leap;

    /* compute leap status for later. 1 = leap, 0 = notleap */
    leap = ((year % 4) == 0) ? 1 : 0;
    if (((year%100) == 0) && ((year%400) != 0)) 
        leap = 0;

    if (month < 1 || month > 12) {
        fprintf(stderr, "error: month must be between 1 and 12\n");
        return ERROR;
    } 
    
    /* check for day range based on month. start with anything but feb */
    if (month != 2  &&  day > days_per_month[month]) {
        fprintf(stderr, "error: day must be between 1 and %d for month %d\n",
                         days_per_month[month], month);
        return ERROR;
    }
    if (month == 2) {
        if (!leap && day > days_per_month[2]) {
            fprintf(stderr, "error: day must be between 1 and %d for month %d\n",
                             days_per_month[month], month);
            return ERROR;
        }
        if (leap && day > days_per_month[13]) {
            fprintf(stderr, "error: day must be between 1 and %d for month %d\n",
                             days_per_month[13], month);
            return ERROR;
        }
    }
    
    if (hour < 0 || hour > 23) {
        fprintf(stderr, "error: hour must be between 0 and 23\n");
        return ERROR;
    } 

    if (min < 0 || min > 59) {
        fprintf(stderr, "error: min must be between 0 and 59\n");
        return ERROR;
    } 

    if (sec < 0 || sec > 59) {
        fprintf(stderr, "error: sec must be between 0 and 59\n");
        return ERROR;
    } 


    return OK;
}


void
convert(int direction, int *days, int *secs, int *year, int *month, int *day,
                       int *hour, int *min, int *sec) {

    int nleapyrs, ndays, nsecs;
    int i, ydays, mdays, leap;

    if (direction == TO_DS) {

        /* compute number of leap years fully past since base_year */
        nleapyrs = (*year - BASE_YEAR) / 4   - 
                   (*year - BASE_YEAR) / 100 + 
                   (*year - BASE_YEAR) / 400;

       /* compute leap status for later. 1 = leap, 0 = notleap */
       leap = ((*year % 4) == 0) ? 1 : 0;
       if (((*year%100) == 0) && ((*year%400) != 0)) 
           leap = 0;

        /* Count up days in this year */
        ndays = 0;
        for (i=1; i<*month; i++) {
            ndays += days_per_month[i];
            if (leap && (i == 2)) ndays += 1;
        }
        
        *secs = *sec + 60 * (*min + 60 * *hour);
        *days = *day - 1 + ndays + 
                 365 * (*year - BASE_YEAR - nleapyrs) + 
                 366 * (nleapyrs);

    } else {
  
        /* if the direction of the loop was reversed, this might work? */
        if (*days < 0) {
            fprintf(stderr, "cannot handle negative days\n");
            return;
        }

        ndays = *days;
        for (i=BASE_YEAR; ; i++) {

            /* Is this a leap year? Gregorian calendar assigns each year evenly
             * divisible by 4 that is not a century year unevenly divisible by 400
             * as a leap-year. (i.e. 1700,1800,1900 are not leap-years, 2000 is)
             */

            /* compute leap status for later. 1 = leap, 0 = notleap */
            leap = ((i % 4) == 0) ? 1 : 0;
            if (((i%100) == 0) && ((i%400) != 0)) 
                leap = 0;

            ydays = 365; 
            if(leap) ydays = 366;
        
            if(ndays >= ydays) 
                ndays -= ydays;
            else {
                *year = i;
                break;
            }
        
        }
        
        /* find month and day */
        for (i=1; i<=12; i++) {
           *month = i;
           mdays = days_per_month[i];
           if(leap && (i == 2)) mdays = 29;
           if(ndays < mdays) break;
           ndays -= mdays;
        }
        
        *day = ndays + 1;
        
        /* Find hour,minute and second */
        nsecs = *secs;
        *hour = nsecs / (60 * 60);
        nsecs -= *hour * (60 * 60);
        *min = nsecs / 60;
        *sec = nsecs - 60 * *min;
        
    }
    
    return;
}

/*
 * <next few lines under version control, do not edit>
 * $URL$
 * $Revision$
 * $Date$
 */
