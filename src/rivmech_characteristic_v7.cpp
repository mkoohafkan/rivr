#include<Rcpp.h>
using namespace Rcpp;

//' @title Froude Number
//' @description Calculate the Froude Number.
//' @param Q Flow rate [\eqn{L^3 T^{-1}}].
//' @param g Gravitational acceleration [\eqn{L T^{-2}}].
//' @param A Flow area [\eqn{L^2}].
//' @param DH Hydraulic depth [\eqn{L}].
//' @return The Froude Number (dimensionless).
//' @details The Froude number is a dimensionless measure of bulk flow 
//'   characteristics that represents the relative importance of inertial 
//'   forces and gravitational forces. For open channel flow, the Froude 
//'   number of open channel flow is defined as \deqn{Fr = \frac{v}{\sqrt{gD_H}}}, 
//'   where \eqn{v = \frac{Q}{A}} is the flow velocity, \eqn{g} is the gravitational 
//'   acceleration and \eqn{D_H} is the hydraulic depth. The Froude number is related
//'   to the energy state of the flow and can be used to identify flows as
//'   either supercritical (\eqn{Fr < 1}) or subcritical (\eqn{Fr > 1}).
//' @export
// [[Rcpp::export]]
double froude(double Q, double g, double A, double DH){
  return((Q/A)/sqrt(DH*g));
}

//' @title Channel geometry
//' @description Compute geometry relations for trapezoidal channels.
//' @param y Flow depth [\eqn{L}].
//' @param B Channel bottom width [\eqn{L}].
//' @param SS Channel sideslope [\eqn{L L^{-1}}]. For a rectangular channel, \code{SS = 0}.
//' @return Named vector:
//'   \item{A}{Flow area [\eqn{L^2}].}
//'   \item{P}{Wetted perimeter [\eqn{L}].}
//'   \item{R}{Hydraulic radius [\eqn{L}].}
//'   \item{dAdy}{Water surface width [\eqn{L}].}
//'   \item{dPdy}{First derivative of wetted perimeter w.r.t. flow depth.}
//'   \item{dRdy}{First derivative of hydraulic radius w.r.t. flow depth.}
//'   \item{DH}{Hydraulic depth [\eqn{L}].}
//'   \item{ybar}{Vertical distance from water surface to centroid of flow area [\eqn{L}].}
//' @details Channel geometry relations are routinely calculated in numerical solutions of steady, 
//'   gradually-varied and unsteady flows. This function is used extensively by internal 
//'   procedures and is made accessible to the user for convenience.
//' @examples
//' channel_geom(1.71, 100, 0) # rectangular channel
//' channel_geom(5.79, 6.1, 1.5) # trapezoidal channel with sideslope 3H:2V
//' @export
// [[Rcpp::export]]
NumericVector channel_geom(double y, double B, double SS){
  // if SS is 0, channel is rectangular
  double A = y*B + SS*y*y; // A
  double P = B + 2*y*sqrt(SS*SS + 1); // P
  double R = A/P; // R
  double dAdy = B + 2*y*SS; // dAdy
  double dPdy = 2*sqrt(SS*SS + 1); // dPdy  
  double dTdy = 2*SS; // dTdy
  double dRdy = dAdy/P - A*dPdy/(P*P); // dRdy
  double DH = A/dAdy;
  double ybar = y*(2*B + dAdy)/(3*(B + dAdy));
  return(NumericVector::create(_["A"] = A, _["P"] = P, _["R"] = R, 
    _["dAdy"] = dAdy, _["dTdy"] = dTdy, _["dPdy"] = dPdy, _["dRdy"] = dRdy, 
    _["DH"] = DH, _["ybar"] = ybar));
}

//' @title Channel conveyance
//' @description Calculate the channel conveyance.
//' @param n Manning's roughness coefficient (dimensionless).
//' @param A Flow area [\eqn{L^2}].
//' @param R Hydraulic radius [\eqn{L}].
//' @param Cm Unit conversion coefficient for Manning's equation. For SI units, Cm = 1.
//' @return The channel conveyance.
//' @details Channel conveyance is routinely calculated in numerical solutions of steady, 
//'   gradually-varied and unsteady flows. This function is used extensively by internal 
//'   procedures and is made accessible to the user for convenience.
//' @export
// [[Rcpp::export]]
double conveyance(double n, double A, double R, double Cm){
  double C = Cm*A*pow(R, 2.0/3.0)/n;
  return(C);
}

double depthfromarea(double A, double yopt, double B, double SS){
  double dy = 9999;
  int maxit = 1000;
  int i = 0;
  double tol = 0.00001;
  while((fabs(dy) > tol)& (i < maxit)){
    NumericVector gp = channel_geom(yopt, B, SS);
    dy = (gp["A"] - A)/gp["dAdy"];
    yopt -= dy;
    i++;
  }
  return(yopt);
}

//' @title Normal depth
//' @description Calculate the normal (equilibrium) depth using Manning's equation.
//' @param So Channel slope [\eqn{L L^{-1}}].
//' @param n Manning's roughness coefficient.
//' @param Q Flow rate [\eqn{L^3 T^{-1}}].
//' @param yopt Initial guess for normal depth [\eqn{L}].
//' @param Cm Unit conversion coefficient for Manning's equation. For SI units, Cm = 1.
//' @param B Channel bottom width [\eqn{L}].
//' @param SS Channel sideslope [\eqn{L L^{-1}}].
//' @return The normal depth \eqn{y_n} [\eqn{L}].
//' @details The normal depth is the equilibrium depth of a channel for a given 
//'   flow rate, channel slope, geometry and roughness.
//'   Manning's equation is used to calculate the equilibrium depth. Manning's 
//'   equation for normal flow is defined as \deqn{Q = \frac{C_m}{n} AR^{2/3}S_0^{1/2}}
//'   where \eqn{Q} is the channel flow, \eqn{S_0} is the channel slope, \eqn{A} is the 
//'   cross-sectional flow area, \eqn{R} is the hydraulic depth and \eqn{C_m} is a conversion factor
//'   based on the unit system used. This function uses a Newton-Raphson root-finding approach 
//'   to calculate the normal depth, i.e.
//'   \eqn{y = y_n} when \deqn{f(y) = \frac{A^{5/3}}{P^{2/3}} - \frac{nQ}{C_mS_0^{1/2}} = 0}.
//' @examples
//' normal_depth(0.001, 0.045, 250, 3, 1.486, 100, 0) # rectangular channel
//' normal_depth(0.0008, 0.013, 126, 5, 1, 6.1, 1.5) # trapezoidal channel with sideslope 3H:2V
//' @export
// [[Rcpp::export]]
double normal_depth(double So, double n, double Q, double yopt, double Cm, 
  double B, double SS){
  // initialize
  double tol = 0.00001;
  int maxit = 1000;
  double dy = 9999;
  int i = 0;
  while((fabs(dy) > tol) & (i < maxit)){
    NumericVector gp = channel_geom(yopt, B, SS) ;
    dy = ( pow(gp["A"], 5.0/3.0)/pow(gp["P"], 2.0/3.0) - n*Q/(Cm*sqrt(So)) ) / 
      ( gp["dAdy"]*(5.0/3.0)*pow(gp["A"]/gp["P"], 2.0/3.0) - 
      (2.0/3.0)*gp["dPdy"]*pow(gp["A"]/gp["P"], 5.0/3.0) );
    yopt = yopt - dy;
    i++;
  }
  return(yopt);  
}

//' @title Critical depth
//' @description Calculate the critical depth.
//' @details The critical depth is the water depth at which a channel 
//'   flow regime will transition from supercritical to subcritical (or vice versa).
//'   Calculation of the critical depth is based on a specific energy formulation, 
//'   i.e. \deqn{E = y + z + \frac{Q^2}{2gB^2y^2}} where \eqn{y} is the flow depth, \eqn{z} is 
//'   the elevation relative to some datum (assumed to be 0), and the last term 
//'   represents kinetic energy. More specifically, the function operates
//'   by finding the point where the derivative of specific energy w.r.t. \eqn{y} is zero, i.e.
//'   \eqn{y = y_c} when \deqn{\frac{dE}{dy} = 1 - \frac{Q^2}{gA^3}\frac{dA}{dy} = 0}.
//' @param Q Flow rate [\eqn{L^3 T^{-1}}].
//' @param yopt Initial guess for normal depth [\eqn{L}].
//' @param g Gravitational acceleration [\eqn{L T^{-2}}].
//' @param B Channel bottom width [\eqn{L}].
//' @param SS Channel sideslope [\eqn{L L^{-1}}].
//' @return The critical depth \eqn{y_c} [\eqn{L}].
//' @examples
//' critical_depth(250, 32.2, 2, 100, 0) # rectangular channel
//' critical_depth(126, 9.81, 1, 6.1, 1.5) # trapezoidal channel with sideslope 3H:2V
//' @export
// [[Rcpp::export]]
double critical_depth(double Q, double yopt, double g, double B, double SS){
  // initialize
  double tol = 0.00001;
  int maxit = 1000;
  double dy = 9999;
  int i = 0;
  while((fabs(dy) > tol) & (i < maxit)){
    NumericVector gp = channel_geom(yopt, B, SS);
    dy = (pow(gp["A"], 3.0)/gp["dAdy"] - Q*Q/g) / 
      (3*pow(gp["A"], 2.0) - pow(gp["A"], 3.0)*gp["dTdy"]/gp["dAdy"]);
    yopt = yopt - dy;
    i++;
  }
  return(yopt);  
}

NumericVector standard_step(double So, double n, double Q, double Cm, double g, 
  double y, double B, double SS, double z, double x, double stepdist){
  // define control section
  NumericVector gp = channel_geom(y, B, SS);
  double E = y + z + 0.5*pow(Q/gp["A"], 2.0)/g;
  double Sf = pow(Q/conveyance(n, gp["A"], gp["R"], Cm), 2.0);
  // define target section
  double znew = z - So*stepdist;
  double xnew = x + stepdist;
  // initial guess for y
  double ynew = y;
  double ylast = 999;
  int i = 0;
  double maxit = 1000;
  double tol = 0.00001;
  while((fabs(ynew - ylast) > tol) & (i < maxit)){
    ylast = ynew ;
    // calculate geometry using yguess	  
    NumericVector geomguess = channel_geom(ylast, B, SS);	
    // calculate Sfguess using geomguess
    double Sfguess = pow(Q/conveyance(n, geomguess["A"], geomguess["R"], Cm), 2.0);
    // calculate average Sf
    double Sfbar = 0.5*(Sfguess + Sf);
    // update ynew
    double j = 0;
	double dy = 9999;
    while((fabs(dy) > tol) & (j < maxit)){
	  NumericVector geomopt = channel_geom(ynew, B, SS);
	  dy = (ynew + 0.5*pow(Q/geomopt["A"], 2.0)/g + znew - E - 
	    Sfbar*fabs(stepdist)) / (1 - geomopt["dAdy"]*Q*Q / 
		(g*pow(geomopt["A"], 3.0)) + 
		fabs(stepdist)*(Sfbar*geomopt["dAdy"]/geomopt["A"] + 
		(2.0/3.0)*Sfbar*geomopt["dRdy"]/geomopt["R"]));
      ynew -= dy;
      j++;
    }
	i++;
  }
  NumericVector geomnew = channel_geom(ynew, B, SS);
  double Sfnew = pow(Q/conveyance(n, geomnew["A"], geomnew["R"], Cm), 2.0);
  double Frnew = froude(Q, g, geomnew["A"], geomnew["DH"]);
  double vnew = Q/geomnew["A"];
  double Enew = znew + ynew + 0.5*vnew*vnew/g;
  return(NumericVector::create(_["x"] = xnew, _["z"] = znew, _["y"] = ynew, 
    _["v"] = vnew, _["A"]=geomnew["A"], _["Sf"] = Sfnew, _["E"] = Enew, _["Fr"] = Frnew));
}

// [[Rcpp::export]]
NumericMatrix loop_step(double So, double n, double Q, double Cm, 
  double g, double y, double B, double SS, double z, double x, 
  double stepdist, double totaldist){
  int numsteps = totaldist/fabs(stepdist) + 1;
  // columns are x, z, y, v, A, Sf, E, Fr
  NumericMatrix pd(numsteps, 8);  
  pd(0, _) = standard_step(So, n, Q, Cm, g, y, B, SS, z, x, 0);
  for(int i = 1; i < numsteps; i++){
    pd(i, _) = standard_step(So, n, Q, Cm, g, pd(i-1, 2), B, SS, pd(i-1, 1), pd(i-1, 0), stepdist);
  }
  return(pd);
}

// [[Rcpp::export]]
List kinematic_wave(double So, double n, double Cm, double g, double B, 
  double SS, int numnodes, NumericVector bc, double ic, double timestep, 
  double spacestep, IntegerVector mpidx, IntegerVector mtidx){
  // initialize
  NumericVector oldflow(numnodes, ic);
  NumericVector olddepth(numnodes, normal_depth(So, n, oldflow[0], 10, Cm, B, SS));
  NumericVector oldarea(numnodes, channel_geom(olddepth[0], B, SS)["A"]); 
  NumericVector oldvelocity(numnodes, oldflow[0]/oldarea[0]);
  // Rcout << "mp = " << mpidx << std::endl;
  // Rcout << "old flow = " << oldflow[0] << std::endl;
  // Rcout << "old depth = " << olddepth[0] << std::endl;
  // monitoring 
  NumericMatrix mpQ(bc.size(), mpidx.size());
  NumericMatrix mpY(bc.size(), mpidx.size());
  NumericMatrix mpV(bc.size(), mpidx.size());
  NumericMatrix mpA(bc.size(), mpidx.size());
  NumericMatrix mtQ(mtidx.size(), numnodes);
  NumericMatrix mtY(mtidx.size(), numnodes);  
  NumericMatrix mtV(mtidx.size(), numnodes);  
  NumericMatrix mtA(mtidx.size(), numnodes);  
  // monitoring references
  IntegerVector mt(bc.size(), -1);
  IntegerVector mp(numnodes, -1);
  
  // NOTE: boundary condition always occupies first column!  
  for(int m = 0; m < mpidx.size(); m++){
    mp[mpidx[m]] = m;
    mpQ(0, m) = oldflow[mpidx[m]];
    mpY(0, m) = olddepth[mpidx[m]];
    mpV(0, m) = oldvelocity[mpidx[m]];	
    mpA(0, m) = oldarea[mpidx[m]];
  }
  // NOTE: initial condition always occupies first row!
  for(int m = 0; m < mtidx.size(); m++){
    mt[mtidx[m]] = m;
  }
  mtQ(0, _) = oldflow;
  mtY(0, _) = olddepth;
  mtV(0, _) = oldvelocity;
  mtA(0, _) = oldarea;
  // loop through time
  for(int k = 1; k < bc.size(); k++){  
    // initialize
    NumericVector newflow(numnodes);
    NumericVector newdepth(numnodes);
    NumericVector newarea(numnodes);
    NumericVector newvelocity(numnodes);	
    // boundary condition
    newflow[0] = bc[k];
    newdepth[0] = normal_depth(So, n, newflow[0], 10, Cm, B, SS);
    newarea[0] = channel_geom(newdepth[0], B, SS)["A"];
    newvelocity[0] = newflow[0]/newarea[0];
    mpQ(k, 0) = newflow[0];
    mpY(k, 0) = newdepth[0];
    mpV(k, 0) = newvelocity[0];
    mpA(k, 0) = newarea[0];
    // loop through space
    for(int i = 1; i < numnodes; i++){
      // compute flow
      newflow[i] = newflow[i-1] - (newarea[i-1] - oldarea[i-1])*spacestep/timestep;
      // compute depth
      double tol = 0.00001;
      int maxit = 1000;
      double dy = 9999;
      int j = 0;
      newdepth[i] = newdepth[i-1];
      while((abs(dy) > tol) & (j < maxit)){
        NumericVector gp = channel_geom(newdepth[i], B, SS) ;
        dy = (gp["A"] - pow(n*newflow[i]/(Cm*sqrt(So)), 3.0/5.0)*pow(gp["P"], 2.0/5.0)) / 
          (gp["dAdy"] - (2.0/5.0)*gp["dPdy"]*pow(n*newflow[i]/(Cm*gp["P"]*sqrt(So)), 3.0/5.0));		  
        newdepth[i] -= dy;
        j++;
      }
      // compute area
      newarea[i] = channel_geom(newdepth[i], B, SS)["A"];
	  newvelocity[i] = newflow[i]/newarea[i];
      // monitor nodes
      if(mp[i] > 0){
        mpQ(k, mp[i]) = newflow[i];
        mpY(k, mp[i]) = newdepth[i];
		mpA(k, mp[i]) = newarea[i];
		mpV(k, mp[i]) = newvelocity[i];
      }
    }
    // monitor times
	if(mt[k] > 0){
      mtQ(mt[k], _) = newflow;
      mtY(mt[k], _) = newdepth;
      mtV(mt[k], _) = newvelocity;
      mtA(mt[k], _) = newarea;
	}
    oldflow = newflow;
    oldarea = newarea;
    olddepth = newdepth;
	oldvelocity = newvelocity;
  }
  List ret;
  ret["monitorpoints"] = mpQ;
  ret["monitortimes"] = mtQ;
  ret["mpdepth"]= mpY;
  ret["mpvelocity"] = mpV;
  ret["mparea"] = mpA;
  ret["mtdepth"] = mtY;
  ret["mtvelocity"] = mtV;
  ret["mtarea"] = mtA;
  return(ret); 
}

double characteristic_velocity(double yn, double vo, double yo, 
  double Sfo, double So, double g, double dt, int curvesign){
  double psi = sqrt(g/yo);
  double phi = vo + curvesign*psi*yo + g*(So - Sfo)*dt;
  double vn = phi - curvesign*psi*yn;
  //Rcout << "yn = " << yn << ", vo = " << vo << ", yo = " << yo << ", Sfo = " << Sfo << ", g = " << g << ", dt = " << dt << ", C = " << curvesign << std::endl;
  return(vn);
}

double characteristic_depth(double vn, double vo, double yo, 
  double Sfo, double So, double g, double dt, int curvesign){
  double psi = sqrt(g/yo);
  double phi = vo + curvesign*psi*yo + g*(So - Sfo)*dt;
  double yn = (phi - vn)/(curvesign*psi);
  //Rcout << "vn = " << vn << ", vo = " << vo << ", yo = " << yo << ", Sfo = " << Sfo << ", g = " << g << ", dt = " << dt << ", C = " << curvesign << std::endl;
  return(yn);
}

double raphson_y(double Qn, double yopt, double B, double SS, double So, double g, 
  double vo, double yo, double Sfo, double dt, int curvesign){
  double tol = 0.00001;
  int maxit = 1000;
  double dy = 9999;
  int i = 0;
  double psi = sqrt(g/yo);
  double phi = vo + curvesign*psi*yo + g*(So - Sfo)*dt;
  while((fabs(dy) > tol) & (i < maxit)){
    NumericVector gp = channel_geom(yopt, B, SS) ;
    dy = (Qn - phi*gp["A"] + curvesign*psi*gp["A"]*yopt) / 
	 (-phi*gp["dAdy"] + curvesign*psi*(gp["A"] + yopt*gp["dAdy"]));
	//Rcout << "dy = " << dy << std::endl;
    yopt -= dy;
    i++;
  }
  return(yopt);
}

// [[Rcpp::export]]
List characteristic_wave(double So, double n, double Cm, double g, double B, 
  double SS, int numnodes, NumericVector bc, NumericVector dc, double ic, 
  double timestep, double spacestep, IntegerVector mpidx, IntegerVector mtidx, 
  std::string btype){
  // initial condition
  NumericVector oldflow(numnodes, ic);
  double idepth;
  if(btype.at(0) == 'y'){
    idepth = bc[0];
  } else {
    idepth =  normal_depth(So, n, oldflow[0], 10, Cm, B, SS);
  }
  NumericVector olddepth(numnodes, idepth);
  NumericVector oldarea(numnodes, channel_geom(olddepth[0], B, SS)["A"]); 
  NumericVector oldconveyance(numnodes, conveyance(n, oldarea[0], 
    channel_geom(olddepth[0], B, SS)["R"], Cm));   
  NumericVector oldsf(numnodes, pow(oldflow[0]/oldconveyance[0], 2.0)); 
  NumericVector oldvelocity(numnodes, oldflow[0]/oldarea[0]);
  NumericVector oldF(numnodes, oldflow[0]*oldvelocity[0] + 
    channel_geom(olddepth[0], B, SS)["ybar"]*oldarea[0]*g);
  NumericVector oldS(numnodes, g*oldarea[0]*(oldsf[0] - So));  
  // initial boundaries
  if(btype.at(0) == 'y'){
    olddepth[0] = bc[0];
  } else {
    oldflow[0] = bc[0];
  }
  if(btype.at(1) == 'y'){
    if(dc[0] < 0){
      olddepth[numnodes-1] = olddepth[numnodes-2];
    } else {
      olddepth[numnodes-1] = dc[0];
    }
  } else {
    if(dc[0] < 0){
      oldflow[numnodes-1] = oldflow[numnodes-2];
    } else {
      oldflow[numnodes-1] = dc[0];
	}
  }  
  // monitoring containers
  NumericMatrix mpQ(bc.size(), mpidx.size());
  NumericMatrix mpY(bc.size(), mpidx.size());
  NumericMatrix mpV(bc.size(), mpidx.size());
  NumericMatrix mpA(bc.size(), mpidx.size());
  NumericMatrix mtQ(mtidx.size(), numnodes);
  NumericMatrix mtY(mtidx.size(), numnodes);  
  NumericMatrix mtV(mtidx.size(), numnodes);  
  NumericMatrix mtA(mtidx.size(), numnodes);  
  // monitoring references
  IntegerVector mt(bc.size(), -1);
  IntegerVector mp(numnodes, -1);
  // NOTE: upstream boundary always occupies first column!
  for(int m = 0; m < mpidx.size(); m++){
    mp[mpidx[m]] = m;
    mpQ(0, m) = oldflow[mpidx[m]];
    mpY(0, m) = olddepth[mpidx[m]];
    mpV(0, m) = oldvelocity[mpidx[m]];	
    mpA(0, m) = oldarea[mpidx[m]];
  }
  // NOTE: initial condition always occupies first row!
  for(int m = 0; m < mtidx.size(); m++){
    mt[mtidx[m]] = m;
  }
  mtQ(0, _) = oldflow;
  mtY(0, _) = olddepth;
  mtV(0, _) = oldvelocity;
  mtA(0, _) = oldarea;
  // loop through time
  for(int k = 1; k < bc.size(); k++){  
    //Rcout << "tstep " << k << std::endl;
    // initialize future variables
    NumericVector newflow(numnodes);
    NumericVector newdepth(numnodes);
    NumericVector newarea(numnodes);
    NumericVector newconveyance(numnodes);   
    NumericVector newsf(numnodes); 
    NumericVector newvelocity(numnodes);
    NumericVector newF(numnodes);
    NumericVector newS(numnodes);
    // upstream boundary (characteristic)
    if(btype.at(0) == 'Q'){ // flow boundary
      newflow[0] = bc[k];
	  if(newflow[0] != 0){
        newdepth[0] = raphson_y(newflow[0], olddepth[0], B, SS, So, 
          g, oldvelocity[1], olddepth[1], oldsf[1], timestep, -1);
	  } else { // solve directly if Q = 0
	    newdepth[0] = characteristic_depth(0, oldvelocity[1], olddepth[1], 
		  oldsf[1], So, g, timestep, -1);
	  }
      //Rcout << "solved depth = " << newdepth[0] << std::endl;
      newarea[0] = channel_geom(newdepth[0], B, SS)["A"];
      newvelocity[0] = newflow[0]/newarea[0];
    } else { // depth boundary
      newdepth[0] = bc[k];
      newvelocity[0] = characteristic_velocity(newdepth[0], oldvelocity[1], 
	    olddepth[1], oldsf[1], So, g, timestep, -1);
      newarea[0] = channel_geom(newdepth[0], B, SS)["A"];
      newflow[0] = newvelocity[0]*newarea[0];
    }
    newconveyance[0] = conveyance(n, newarea[0], 
      channel_geom(newdepth[0], B, SS)["R"], Cm);
    newsf[0] = pow(newflow[0]/newconveyance[0], 2.0);
    newF[0] = newflow[0]*newvelocity[0] + 
	  channel_geom(newdepth[0], B, SS)["ybar"]*newarea[0]*g;
    newS[0] = g*newarea[0]*(newsf[0] - So);
    // monitor node
    mpQ(k, 0) = newflow[0];
    mpY(k, 0) = newdepth[0];
    mpV(k, 0) = newvelocity[0];
    mpA(k, 0) = newarea[0];
    // loop through space: predictor
    NumericVector flowstar(numnodes);
    NumericVector areastar(numnodes);
    NumericVector depthstar(numnodes);
    NumericVector conveyancestar(numnodes);
    NumericVector sfstar(numnodes);	
    NumericVector velocitystar(numnodes);
    NumericVector Fstar(numnodes);
    NumericVector Sstar(numnodes);
	NumericVector areastargrad(numnodes);
	NumericVector flowstargrad(numnodes);	
	// loop through space: predictor
    for(int i = 1; i < numnodes; i++){	
	  flowstargrad[i] = - timestep*(oldF[i] - oldF[i-1])/spacestep - 
        timestep*oldS[i];
      flowstar[i] = oldflow[i] + flowstargrad[i];
	  areastargrad[i] = - timestep*(oldflow[i] - oldflow[i-1])/spacestep;
      areastar[i] = oldarea[i] + areastargrad[i];
      depthstar[i] = depthfromarea(areastar[i], olddepth[i], B, SS);    
      conveyancestar[i] = conveyance(n, areastar[i], channel_geom(depthstar[i], 
        B, SS)["R"], Cm);
      sfstar[i] = pow(flowstar[i]/conveyancestar[i], 2.0);
      velocitystar[i] = flowstar[i]/areastar[i];
      Fstar[i] = flowstar[i]*velocitystar[i] + 
	    g*channel_geom(depthstar[i], B, SS)["ybar"]*areastar[i];
      Sstar[i] = g*areastar[i]*(sfstar[i] - So);
    }
    // loop through space: corrector
	NumericVector areadstargrad(numnodes);
	NumericVector flowdstargrad(numnodes);	
    for(int i = 1; i < numnodes - 1; i++){	
      flowdstargrad[i] =  - timestep*(Fstar[i+1] - Fstar[i])/spacestep - timestep*Sstar[i];
      areadstargrad[i] = - timestep*(flowstar[i+1] - flowstar[i])/spacestep;
    }
    // loop through space: final
    for(int i = 1; i < numnodes - 1; i++){
      newflow[i] = oldflow[i] + 0.5*(flowstargrad[i] + flowdstargrad[i]);
      newarea[i] = oldarea[i] + 0.5*(areastargrad[i] + areadstargrad[i]);
      newvelocity[i] = newflow[i]/newarea[i];
      newdepth[i] = depthfromarea(newarea[i], olddepth[i], B, SS);
	  //Rcout << "A* = " << areastar[i] <<  ", A** = " << areadstar[i] << ", area = "<< newarea[i] << ", y = " << newdepth[i] << std::endl;
      newconveyance[i] = conveyance(n, newarea[i], 
        channel_geom(newdepth[i], B, SS)["R"], Cm);
      newsf[i] = pow(newflow[i]/newconveyance[i], 2);
      newF[i] = newflow[i]*newvelocity[i] + 
	    g*channel_geom(newdepth[i], B, SS)["ybar"]*newarea[i];
      newS[i] = g*newarea[i]*(newsf[i] - So);
      // monitor nodes
      if(mp[i] >= 0){
        mpQ(k, mp[i]) = newflow[i];
        mpY(k, mp[i]) = newdepth[i];
		mpA(k, mp[i]) = newarea[i];
		mpV(k, mp[i]) = newvelocity[i];
      }
    }
    // downstream boundary (characteristic)
    if(btype.at(1) == 'Q'){ // flow boundary
      if(dc[k] < 0){ // use flow of adjacent node if dc[k] < 0 (zero gradient)
        newflow[numnodes-1] = newflow[numnodes-2];
      } else {
        newflow[numnodes-1] = dc[k];
      }
	  if(newflow[numnodes-1] != 0){
	    //Rcout << "newflow[numnodes-1] = " << newflow[numnodes-1] << std::endl;
        newdepth[numnodes-1] = raphson_y(newflow[numnodes-1],olddepth[numnodes-1], 
          B, SS, So, g, oldvelocity[numnodes-2], olddepth[numnodes-2], 
          oldsf[numnodes-2], timestep, 1);
	  } else { // solve directly if flow is zero
	    newdepth[numnodes-1] = characteristic_depth(0, oldvelocity[numnodes-2], 
		  olddepth[numnodes-2], oldsf[numnodes-2], So, g, timestep, 1);
	  }
	  //Rcout << "v(" << k << ", 1end) = " << newvelocity[numnodes-1] << "; y(" << k << ", end) = " << newdepth[numnodes-1] << std::endl;
      newarea[numnodes-1] = channel_geom(newdepth[numnodes-1], B, SS)["A"];
      newvelocity[numnodes-1] = newflow[numnodes-1]/newarea[numnodes-1];
    } else {
      if(dc[k] < 0){
        newdepth[numnodes-1] = newdepth[numnodes-2];
      } else {
        newdepth[numnodes-1] = dc[k];
      }
      newvelocity[numnodes-1] = characteristic_velocity(newdepth[numnodes-1], 
        oldvelocity[numnodes-2], olddepth[numnodes-2], oldsf[numnodes-2], 
        So, g, timestep, 1);
      newarea[numnodes-1] = channel_geom(newdepth[numnodes-1], B, SS)["A"];
      newflow[numnodes-1] = newvelocity[numnodes-1]*newarea[numnodes-1];
	}
    newconveyance[numnodes-1] = conveyance(n, newarea[numnodes-1], 
        channel_geom(newdepth[numnodes-1], B, SS)["R"], Cm);		
    newsf[numnodes-1] = pow(newflow[numnodes-1]/newconveyance[numnodes-1], 2);
    newF[numnodes-1] = newflow[numnodes-1]*newvelocity[numnodes-1] + 
      g*channel_geom(newdepth[numnodes-1], B, SS)["ybar"]*newarea[numnodes-1];
    newS[numnodes-1] = g*newarea[numnodes-1]*(newsf[numnodes-1] - So);
    //Rcout << "v(" << k << ", 1end) = " << newvelocity[numnodes-1] << "; y(" << k << ", end) = " << newdepth[numnodes-1] << std::endl;
    // monitor nodes
    if(mp[numnodes-1] >= 0){
      mpQ(k, mp[numnodes-1]) = newflow[numnodes-1];
      mpY(k, mp[numnodes-1]) = newdepth[numnodes-1];
      mpV(k, mp[numnodes-1]) = newvelocity[numnodes-1];
      mpA(k, mp[numnodes-1]) = newarea[numnodes-1];  
    }
    // monitor times
    if(mt[k] >= 0){
      mtQ(mt[k], _) = newflow;
      mtY(mt[k], _) = newdepth;
      mtV(mt[k], _) = newvelocity;
      mtA(mt[k], _) = newarea;
    }
    oldflow = newflow;
    oldarea = newarea;
    olddepth = newdepth;
    oldconveyance = newconveyance;   
    oldsf = newsf; 
    oldvelocity = newvelocity;
    oldF = newF;
    oldS = newS;  
  }
  //Rcout << "last upstream flow = " << newflow[0] << std::endl;
  //Rcout << "last upstream depth = " << newdepth[0] << std::endl;
  //Rcout << "last upstream velocity = " << newvelocity[0] << std::endl;
  List ret;
  ret["monitorpoints"] = mpQ;
  ret["monitortimes"] = mtQ;
  ret["mpdepth"]= mpY;
  ret["mtdepth"] = mtY;  
  ret["mpvelocity"] = mpV;
  ret["mtvelocity"] = mtV;
  ret["mparea"] = mpA;
  ret["mtarea"] = mtA;
  return(ret); 
}

// [[Rcpp::export]]
List diffusive_wave(double So, double n, double Cm, double g, double B, 
  double SS, int numnodes, NumericVector bc, NumericVector dc, double ic, 
  double timestep, double spacestep, IntegerVector mpidx, IntegerVector mtidx, 
  std::string btype){
  // initial condition
  NumericVector oldflow(numnodes, ic);
  double idepth;
  if(btype.at(0) == 'y'){
    idepth = bc[0];
  } else {
    idepth =  normal_depth(So, n, oldflow[0], 10, Cm, B, SS);
  }
  NumericVector olddepth(numnodes, idepth);
  NumericVector oldarea(numnodes, channel_geom(olddepth[0], B, SS)["A"]); 
  NumericVector oldconveyance(numnodes, conveyance(n, oldarea[0], 
    channel_geom(olddepth[0], B, SS)["R"], Cm));   
  NumericVector oldsf(numnodes, pow(oldflow[0]/oldconveyance[0], 2.0)); 
  NumericVector oldvelocity(numnodes, oldflow[0]/oldarea[0]);
  NumericVector oldF(numnodes, oldflow[0]*oldvelocity[0] + 
    channel_geom(olddepth[0], B, SS)["ybar"]*oldarea[0]*g);
  NumericVector oldS(numnodes, g*oldarea[0]*(oldsf[0] - So));  
  // monitoring containers
  NumericMatrix mpQ(bc.size(), mpidx.size());
  NumericMatrix mpY(bc.size(), mpidx.size());
  NumericMatrix mpV(bc.size(), mpidx.size());
  NumericMatrix mpA(bc.size(), mpidx.size());
  NumericMatrix mtQ(mtidx.size(), numnodes);
  NumericMatrix mtY(mtidx.size(), numnodes);  
  NumericMatrix mtV(mtidx.size(), numnodes);  
  NumericMatrix mtA(mtidx.size(), numnodes);  
  // monitoring references
  IntegerVector mt(bc.size(), -1);
  IntegerVector mp(numnodes, -1);
  // NOTE: upstream boundary always occupies first column!
  for(int m = 0; m < mpidx.size(); m++){
    mp[mpidx[m]] = m;
    mpQ(0, m) = oldflow[mpidx[m]];
    mpY(0, m) = olddepth[mpidx[m]];
    mpV(0, m) = oldvelocity[mpidx[m]];	
    mpA(0, m) = oldarea[mpidx[m]];
  }
  // NOTE: initial condition always occupies first row!
  for(int m = 0; m < mtidx.size(); m++){
    mt[mtidx[m]] = m;
  }
  mtQ(0, _) = oldflow;
  mtY(0, _) = olddepth;
  mtV(0, _) = oldvelocity;
  mtA(0, _) = oldarea;
  // loop through time
  for(int k = 1; k < bc.size(); k++){  
    //Rcout << "tstep " << k << std::endl;
    // initialize future variables
    NumericVector newflow(numnodes);
    NumericVector newdepth(numnodes);
    NumericVector newarea(numnodes);
    NumericVector newconveyance(numnodes);   
    NumericVector newsf(numnodes); 
    NumericVector newvelocity(numnodes);
    NumericVector newF(numnodes);
    NumericVector newS(numnodes);
    // upstream boundary (characteristic)
    if(btype.at(0) == 'Q'){ // flow boundary
      newflow[0] = bc[k];
	  if(newflow[0] != 0){
        newdepth[0] = raphson_y(newflow[0], olddepth[0], B, SS, So, 
          g, oldvelocity[1], olddepth[1], oldsf[1], timestep, -1);
	  } else { // solve directly if Q = 0
	    newdepth[0] = characteristic_depth(0, oldvelocity[1], olddepth[1], 
		  oldsf[1], So, g, timestep, -1);
	  }
      //Rcout << "solved depth = " << newdepth[0] << std::endl;
      newarea[0] = channel_geom(newdepth[0], B, SS)["A"];
      newvelocity[0] = newflow[0]/newarea[0];
    } else { // depth boundary
      newdepth[0] = bc[k];
      newvelocity[0] = characteristic_velocity(newdepth[0], oldvelocity[1], 
	    olddepth[1], oldsf[1], So, g, timestep, -1);
      newarea[0] = channel_geom(newdepth[0], B, SS)["A"];
      newflow[0] = newvelocity[0]*newarea[0];
    }
    newconveyance[0] = conveyance(n, newarea[0], 
      channel_geom(newdepth[0], B, SS)["R"], Cm);
    newsf[0] = pow(newflow[0]/newconveyance[0], 2.0);
    newF[0] = newflow[0]*newvelocity[0] + 
	  channel_geom(newdepth[0], B, SS)["ybar"]*newarea[0]*g;
    newS[0] = g*newarea[0]*(newsf[0] - So);
    // downstream boundary (characteristic)
    if(btype.at(1) == 'Q'){ // flow boundary
      if(dc[k] < 0){ // use flow of adjacent node if dc[k] < 0 (zero gradient)
        newflow[numnodes-1] = newflow[numnodes-2];
      } else {
        newflow[numnodes-1] = dc[k];
      }
	  if(newflow[numnodes-1] != 0){
	    //Rcout << "newflow[numnodes-1] = " << newflow[numnodes-1] << std::endl;
        newdepth[numnodes-1] = raphson_y(newflow[numnodes-1],olddepth[numnodes-1], 
          B, SS, So, g, oldvelocity[numnodes-2], olddepth[numnodes-2], 
          oldsf[numnodes-2], timestep, 1);
	  } else { // solve directly if flow is zero
	    newdepth[numnodes-1] = characteristic_depth(0, oldvelocity[numnodes-2], 
		  olddepth[numnodes-2], oldsf[numnodes-2], So, g, timestep, 1);
	  }
	  //Rcout << "v(" << k << ", 1end) = " << newvelocity[numnodes-1] << "; y(" << k << ", end) = " << newdepth[numnodes-1] << std::endl;
      newarea[numnodes-1] = channel_geom(newdepth[numnodes-1], B, SS)["A"];
      newvelocity[numnodes-1] = newflow[numnodes-1]/newarea[numnodes-1];
    } else {
      if(dc[k] < 0){
        newdepth[numnodes-1] = newdepth[numnodes-2];
      } else {
        newdepth[numnodes-1] = dc[k];
      }
      newvelocity[numnodes-1] = characteristic_velocity(newdepth[numnodes-1], 
        oldvelocity[numnodes-2], olddepth[numnodes-2], oldsf[numnodes-2], 
        So, g, timestep, 1);
      newarea[numnodes-1] = channel_geom(newdepth[numnodes-1], B, SS)["A"];
      newflow[numnodes-1] = newvelocity[numnodes-1]*newarea[numnodes-1];
	}
    newconveyance[numnodes-1] = conveyance(n, newarea[numnodes-1], 
        channel_geom(newdepth[numnodes-1], B, SS)["R"], Cm);		
    newsf[numnodes-1] = pow(newflow[numnodes-1]/newconveyance[numnodes-1], 2);
    newF[numnodes-1] = newflow[numnodes-1]*newvelocity[numnodes-1] + 
      g*channel_geom(newdepth[numnodes-1], B, SS)["ybar"]*newarea[numnodes-1];
    newS[numnodes-1] = g*newarea[numnodes-1]*(newsf[numnodes-1] - So);
    // monitor node
    mpQ(k, 0) = newflow[0];
    mpY(k, 0) = newdepth[0];
    mpV(k, 0) = newvelocity[0];
    mpA(k, 0) = newarea[0];	
    // loop through space
    for(int i = 1; i < numnodes - 1; i++){	  
	  double Astar = 0.5*(oldarea[i+1] + oldarea[i-1]);
      double flowstar = 0.5*(oldflow[i+1] + oldflow[i-1]);
	  double Sstar = 0.5*(oldS[i+1] + oldS[i-1]);
      newarea[i] = Astar - 0.5*timestep*(oldflow[i+1] - oldflow[i-1])/spacestep;
	  newflow[i] = flowstar - 0.5*timestep*(oldF[i+1] - oldF[i-1])/spacestep - timestep*Sstar;
	  newvelocity[i] = newflow[i]/newarea[i];
	  newdepth[i] = depthfromarea(newarea[i], olddepth[i], B, SS);
	  newconveyance[i] = conveyance(n, newarea[i], 
        channel_geom(newdepth[i], B, SS)["R"], Cm);
	  newsf[i] = pow(newflow[i]/newconveyance[i], 2);
      newF[i] = newflow[i]*newvelocity[i] + 
	    g*channel_geom(newdepth[i], B, SS)["ybar"]*newarea[i];
	  newS[i] = g*newarea[i]*(newsf[i] - So);
      // monitor nodes
      if(mp[i] >= 0){
        mpQ(k, mp[i]) = newflow[i];
        mpY(k, mp[i]) = newdepth[i];
		mpA(k, mp[i]) = newarea[i];
		mpV(k, mp[i]) = newvelocity[i];
      }
    }
    // monitor nodes
    if(mp[numnodes-1] >= 0){
      mpQ(k, mp[numnodes-1]) = newflow[numnodes-1];
      mpY(k, mp[numnodes-1]) = newdepth[numnodes-1];
      mpV(k, mp[numnodes-1]) = newvelocity[numnodes-1];
      mpA(k, mp[numnodes-1]) = newarea[numnodes-1];  
    }
    // monitor times
    if(mt[k] >= 0){
      mtQ(mt[k], _) = newflow;
      mtY(mt[k], _) = newdepth;
      mtV(mt[k], _) = newvelocity;
      mtA(mt[k], _) = newarea;
    }
    oldflow = newflow;
    oldarea = newarea;
    olddepth = newdepth;
    oldconveyance = newconveyance;   
    oldsf = newsf; 
    oldvelocity = newvelocity;
    oldF = newF;
    oldS = newS;  
  }
  //Rcout << "last upstream flow = " << newflow[0] << std::endl;
  //Rcout << "last upstream depth = " << newdepth[0] << std::endl;
  //Rcout << "last upstream velocity = " << newvelocity[0] << std::endl;
  List ret;
  ret["monitorpoints"] = mpQ;
  ret["monitortimes"] = mtQ;
  ret["mpdepth"]= mpY;
  ret["mtdepth"] = mtY;  
  ret["mpvelocity"] = mpV;
  ret["mtvelocity"] = mtV;
  ret["mparea"] = mpA;
  ret["mtarea"] = mtA;
  return(ret); 
}
