/* NUMBER THEROY */

// Extended Euclides  a*s + b*t = gcd(a,b)
// finds 'd' and any 's' and 't' 

lli d,s,t;
void gcd(lli a, lli b){
   if(!b){ d=abs(a); s=a>0?1:-1; t=0; }
   else{ gcd(b,a%b); swap(s,t); t=t-(a/b)*s; }
}


// lucas(n,k,p): Binomial Coefficient for  n,k < 10^18  modulo prime  p < 10^7
// Complexity: equal of the complexity of C(n,k,p) 
// C(n,k,p): must return the binomial coefficient nCk mod p, with n,k<p with p prime

lli lucas(lli n, lli k, lli p){
	lli s=1;
	while(n || k){ s=s*C(n%p,k%p,p)%p; n/=p; k/=p;}
	return s;
}


// Miller Rabin Primality Test

#define lli long long

inline lli mul(lli a, lli b, lli M ){ return (__int128(a)*b)%M; }

lli exp(lli b, lli e, lli M){lli s=1; for(;e;e>>=1){ if(e&1) s=mul(s,b,M); b=mul(b,b,M); } return s; }

int prime(lli p, int it=10){
   if(p<2) return 0;
   if(p!=2 && p%2==0) return 0;
   lli s=p-1;
   while(s%2==0){ s/=2; }
   for(int i=0;i<it;i++){
      lli a=rand()%(p-1)+1,tmp=s;
      lli m=exp(a,tmp,p);
      while(tmp!=p-1 && m!=1 && m!=p-1){ m=mul(m,m,p); tmp*=2; }
      if(m!=p-1 && tmp%2==0) return 0;
   }
   return 1;
}

// Pollard Rho using Floyd's Cycle Finding
// Complexity O(n^1/4)
// Returns any factor of 'n'
// if 'n' is prime never returns

lli gcd(lli a,lli b){return b? gcd(b,a%b):abs(a); }

lli rho(lli n){
   lli a=rand(),b=rand(),x=2,y=2,d=1;
   #define f(x) (mul(x,x+a,n)+b)%n
   while(d==1){ x=f(x); y=f(f(y)); d=gcd(x-y,n); }
   return n==d? rho(n):d;
}



/* GEOMETRY */


// Sweep Line for Closest Pair
// Complexity O( nlog(n) ) , n=|pts|
// Finds minimum distance between any pair in pts

#define lli double 
#define pt pair<lli,lli>
#define x first 
#define y second

struct cmp{ bool operator()(pt a, pt b){ return a.y<b.y; } };

lli dist(vector<pt> &pts){
   sort(pts.begin(),pts.end()); set<pt,cmp> sw;
   lli d=10e17; int w=0;
   for(auto p:pts){
      while( abs(p.x-pts[w].x)>d ) sw.erase(pts[w++]); 
      auto it=sw.lower_bound(pt(0,p.y-d));
      auto fin=sw.upper_bound(pt(0,p.y+d));
      for(;it!=fin;++it) d=min(d,hypot(it->x-p.x,it->y-p.y)); 
      sw.insert(p); 
   }
   return d;
}

// Convex Hull 
// Complexity O( nlog(n) ) , n=|pts|
// Return a subset of 'pts' with the points in the
// convex hull sorted in counter clockwise order

lli ccw(pt a, pt b, pt p){ return (b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x); }

vector<pt> hull(vector<pt> pts){
   sort(pts.begin(),pts.end());
   vector<pt> up,dn;
   for(int i=0;i<pts.size();i++){
      while( dn.size()>=2 && ccw(dn[dn.size()-2],dn[dn.size()-1],pts[i])<=0 ) dn.pop_back();
      dn.push_back(pts[i]);
      while( up.size()>=2 && ccw(up[up.size()-2],up[up.size()-2],pts[i])>=0 ) up.pop_back();
      up.push_back(pts[i]);
   }
   dn.insert(dn.end(),up.rbegin()+1,up.rend()-1);
   return dn;
}


/* OTHERS */

// Convert dates to ints

int date_to_int(int d, int m, int y){  
  return 
    1461 * (y + 4800 + (m - 14) / 12) / 4 +
    367 * (m - 2 - (m - 14) / 12 * 12) / 12 - 
    3 * ((y + 4900 + (m - 14) / 12) / 100) / 4 + 
    d - 32075;
}

void int_to_date(int jd, int &d, int &m, int &y){
  int x, n, i, j;
  x = jd + 68569;
  n = 4 * x / 146097;
  x -= (146097 * n + 3) / 4;
  i = (4000 * (x + 1)) / 1461001;
  x -= 1461 * i / 4 - 31;
  j = 80 * x / 2447;
  d = x - 2447 * j / 80;
  x = j / 11;
  m = j + 2 - 12 * x;
  y = 100 * (n - 49) + i + x;
}