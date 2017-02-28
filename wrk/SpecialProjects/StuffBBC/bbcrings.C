int bbcrings(int i)
{

  if ( i == 8 ||
       i == 11||
       i == 14||
       i == 19||
       i == 22||
       i == 26||
       i == 40||
       i == 43||
       i == 46||
       i == 51||
       i == 54||
       i == 58 ) return 0; // inner layer

  if ( i == 7 ||
       i == 16||
       i == 25||
       i == 39||
       i == 48||
       i == 57 ) return 1; // inner middle layer

  if ( i == 4 ||
       i == 6 ||
       i == 10||
       i == 13||
       i == 18||
       i == 21||
       i == 24||
       i == 29||
       i == 36||
       i == 38||
       i == 42||
       i == 45||
       i == 45||
       i == 50||
       i == 53||
       i == 56||
       i == 61 ) return 2; // middle layer

  if ( i == 3 ||
       i == 15||
       i == 28||
       i == 35||
       i == 47||
       i == 60 ) return 3; //outer middle layer

  if ( i == 0 ||
       i == 1 ||
       i == 2 ||
       i == 5 ||
       i == 9 ||
       i == 12||
       i == 17||
       i == 20||
       i == 23||
       i == 27||
       i == 30||
       i == 31||
       i == 32||
       i == 33||
       i == 34||
       i == 37||
       i == 41||
       i == 44||
       i == 49||
       i == 52||
       i == 55||
       i == 59||
       i == 62||
       i == 63 ) return 4; // outer layer

  return -1;

}
