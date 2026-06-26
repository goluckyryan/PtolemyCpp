/*
 *    real*8 dtime
 *    time = dtime()
 */

#include <sys/time.h>
#include <sys/resource.h>
#include <cstdio>

double dtime_( void )
{
     struct rusage myusage;
     double userTime, sysTime;

     if ( getrusage(RUSAGE_SELF, &myusage) < 0 ) {
       perror( "getrusage error" );
       return 0.0;
     }
     userTime = static_cast<double>(myusage.ru_utime.tv_sec) +
	myusage.ru_utime.tv_usec / 1000000.0;
     sysTime    = static_cast<double>(myusage.ru_stime.tv_sec) +
	myusage.ru_stime.tv_usec / 1000000.0;
     return userTime + sysTime;
}
