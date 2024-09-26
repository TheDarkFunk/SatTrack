package SatTrak
  model Satellite
    constant Real pi = Modelica.Constants.pi;
    constant Real d2r = Modelica.Constants.D2R;
    constant Real k2 = 2.2005645e4 "0.5J2Re^2 (km2)";
    constant Real mu = 398600.4;
    parameter Real M0 "Mean anomaly at Epoch (deg)";
    parameter Real N0 "Mean motion at Epoch (rev/d)";
    parameter Real eccn "Eccentricity";
    parameter Real Ndot2 "1st Der Mean Motion over 2 (rev/d2)";
    parameter Real Nddot6 "2nd Der Mean Motion over 6 (rev/d3)";
    parameter Real raan0 "Right Ascension Ascending Node at Epoch (deg)";
    parameter Real argper0 "Argument of Perigee at Epoch (deg)";
    parameter Real incl "Inclination angle (deg)";
    parameter Real tstart "Time from ref epoch to start of sim (sec)";
    Real M "Mean Anomaly (deg)";
    Real N "Mean Motion (rev/d)";
    Real E "Eccentric Anomaly (deg)";
    Real theta "true anomaly (deg)";
    Real a "Semi-major axis (km)";
    Real a0 = (398600.4 * 86400 ^ 2 / 4 / pi ^ 2 / N0 ^ 2) ^ (1 / 3) "Semi-major axis at Epoch";
    Real raan "Right Ascension Ascending Node (deg)";
    Real argper "Argument of Perigee (deg)";
    Real r "Satellite radial distance (km)";
    Real rotang[3] "TODO: fill in the angles here";
    Real p_sat_pf[3] "Satellite posn in pf coords (km)";
    Real v_sat_pf[3] "Satellite vel in pf coords (km/s)";
  initial equation
    M = M0 + N0 * (360. * tstart) / 86400. + 360. * Ndot2 * (tstart / 86400.) ^ 2 + 360. * Nddot6 * (tstart / 86400.) ^ 3;
    N = N0 + 2 * Ndot2 * (time / 86400. ^ 2) + 3 * Nddot6 * (time ^ 2 / 86400. ^ 3);
    a = a0;
    raan = raan0 + N0 * 360. / 86400. * (3. * k2 * cos(incl * d2r) / (a ^ 2 * (1. - eccn ^ 2) ^ 2)) * tstart;
    argper = argper0 + N0 * 360. / 86400. * (-3. * k2 * (5. * cos(incl * d2r) ^ 2 - 1.) / (a ^ 2 * (1. - eccn ^ 2) ^ 2)) * tstart;
    incl = incl;
  equation
    M * d2r = E * d2r - eccn * sin(E * d2r);
    der(M) = N0 * (360. / 86400.) + 2 * 360. * Ndot2 * (time / 86400. ^ 2) + 3 * 360. * Nddot6 * (time ^ 2 / 86400. ^ 3);
    tan(theta * d2r / 2.) = sqrt((1. + eccn) / (1. - eccn)) * tan(E * d2r / 2.);
    r = a * (1. - eccn ^ 2) / (1. + eccn * cos(theta * d2r));
    N = 86400 / (2. * pi) * sqrt(mu / a ^ 3);
    der(N) = 2 * 360. * Ndot2 * (1 / 86400. ^ 2) + 6 * 360. * Nddot6 * (time / 86400. ^ 3);
    p_sat_pf[1] = r * cos(theta * d2r);
    p_sat_pf[2] = r * sin(theta * d2r);
    p_sat_pf[3] = 0.;
    der(p_sat_pf[1]) = v_sat_pf[1];
    der(p_sat_pf[2]) = v_sat_pf[2];
    der(p_sat_pf[3]) = v_sat_pf[3];
    der(raan) = N * 360. / 86400. * (3. * k2 * cos(incl * d2r) / (a ^ 2 * (1 - eccn ^ 2) ^ 2));
    der(argper) = N * 360. / 86400. * (-3. * k2 * (5. * cos(incl * d2r) ^ 2 - 1) / (a ^ 2 * (1 - eccn ^ 2) ^ 2));
    rotang[1] = -argper / d2r;
    rotang[2] = -incl / d2r;
    rotang[3] = -raan / d2r;
  end Satellite;

  model SatTest
    //GPS BIIF-2  (PRN 01)
    //1 37753U 11036A   19066.16055598 -.00000005  00000-0  00000+0 0  9992
    //2 37753  55.8450  75.7803 0083249  37.8271 138.9022  2.00565768 55934
    //instantiate the Satellite model with TLE parameters above
    SatTrak.Satellite SatTestP1(M0, N0, eccn, Ndot2, raan0, Nddot6, argper0, incl, tstart);
    //instantiate the Station model with ARO parameters
    SatTrak.Station SatTestP2(stn_long, stn_lat, stn_elev);
    //Trajectory Parameters
    constant Real d2r = Modelica.Constants.D2R;
    parameter Real hour0 "dayfraction at starting epoch";
    parameter Real day "Number of days from J2000 to midnight of epoch day";
    Real M "Mean Anomaly (deg)";
    Real N "Mean Motion (rev/d)";
    Real E "Eccentric Anomaly (deg)";
    Real theta "true anomaly (deg)";
    Real incl "inclination (deg)";
    Real a "Semi-major axis (km)";
    Real raan "Right Ascension Ascending Node (deg)";
    Real argper "Argument of Perigee (deg)";
    Real r "Satellite radial distance (km)";
    Real p_sat_pf[3] "Satellite posn in pf coords (km)";
    Real v_sat_pf[3] "Satellite vel in pf coords (km/s)";
    Real rotang[3] "Rotation Angles";
    //P2 parameters
    Real p_sat_ECI[3] "Satellite position in ECI coordinates (km)";
    Real v_sat_ECI[3] "Satellite velocity in ECI coordinates (km/s)";
    Real GMST "GMST angle (deg)";
    Real p_sat_ECF[3] "Position vector in ECF coordinates (km)";
    Real v_sat_ECF[3] "Relative Velocity vector in ECF coordinates (km/s)";
    Real p_sat_topo[3] "position of satellite relative to station (km)";
    Real v_sat_topo[3] "velocity of satellite relative to sration (km/s)";
    Real p_stn_topo[3] "Station coordinates in topocentric (km)";
    Real p_stn_ECF[3] "Station coordinates in ECF (km)";
    Real hours "day fraction at simulation time";
    Real TM[3, 3] "Transform matrix from ECF to topo";
    Real azimuth;
    Real elevation;
    Real azdot;
    Real eldot;
    Real Elmin = 9 "Minimum Elevation For Station (deg)";
    Real Elmax = 89 "Maximum Elevation for Station (deg)";
    Real AOS "Acqu of signal Time (epoch s)";
    Real LOS "Loss of signal Time (epoch s)";
    Boolean InsideFOR "Satellite is within FOR";
    Boolean OutsideFOR "Satellite is outside FOR (Logic Inverse of InsideFOR)";
  equation
//TEST the Satelite model (simulating trajectory)
    M = mod(SatTestP1.M, 360);
    N = SatTestP1.N;
    E = SatTestP1.E;
    theta = SatTestP1.theta;
    a = SatTestP1.a;
    incl = SatTestP1.incl;
    raan = SatTestP1.raan;
    argper = SatTestP1.argper;
    r = SatTestP1.r;
    p_sat_pf = SatTestP1.p_sat_pf;
    v_sat_pf = SatTestP1.v_sat_pf;
    rotang = SatTestP1.rotang;
//Beginning of P2 testing
//TEST perifocal to earth centered inertial conversion function
    (p_sat_ECI, v_sat_ECI) = Sat_ECI(ang = {-1 * d2r * argper, -1 * d2r * incl, -1 * d2r * raan}, p_pf = p_sat_pf, v_pf = v_sat_pf);
//TEST GMST angle from epoch date
//recall that from TLE 19066.16055598 gives year (2019) and day of year (66.16055598)
// theta_d requires day input to be days elapsed since J2000
// so Day = 66+2458484 (2019) - 2451544 (J2000) and Hour = 0.16055598
    hours = hour0 + time * (1 / 86400) "start dayfrac + simulation time (seconds converted to dayfrac)";
    GMST = theta_d(days = day, hours = hours);
//TEST satellite ECI to ECF coordinates
    (p_sat_ECF, v_sat_ECF) = sat_ECF(ang = GMST, p_ECI = p_sat_ECI, v_ECI = v_sat_ECI);
//TEST range function, trajectory of satellite relative to station
    p_stn_ECF = SatTestP2.p_stn_ECF "get station ECF position from station model";
    p_stn_topo = SatTestP2.p_stn_top "station topo position from Station model";
    TM = SatTestP2.TM "get ECF to topo transformation from station model";
    (p_sat_topo, v_sat_topo) = range_ECF2topo(p_stn_ECF, p_sat_ECF, v_sat_ECF, TM);
//TEST look angle function
    (azimuth, elevation, azdot, eldot) = range_topo2look_angles(p_sat_topo = p_sat_topo, v_sat_topo = v_sat_topo);
//TEST Visibility
    InsideFOR = Elmin <= elevation and elevation <= Elmax;
    OutsideFOR = elevation > Elmax or elevation < Elmin;
    if initial() then
      if InsideFOR then
        AOS = -10 "cases 2 and 3 starting inside FOR";
      end if;
    end if;
//edge is always from False to True
    when edge(InsideFOR) then
      AOS = time "Case 1 when entering FOR";
    end when;
    when edge(OutsideFOR) then
      LOS = time "Case 3 when exiting FOR";
    end when;
  end SatTest;

  model Station
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    constant Real Re = 6378.137 "Earth radius (km)";
    constant Real f = 1 / 298.257223563 "Earth ref ellipsoid flattening";
    constant Real ecc_e_sq = 2 * f - f ^ 2 "Square of earth eccentricity";
    constant Real d2r = Modelica.Constants.D2R;
    parameter Real stn_long "Station longitude (degE)";
    parameter Real stn_lat "Station latitude (degN)";
    parameter Real stn_elev "Station elevation (m)";
    Real p_stn_ECF[3] "Station coordinates in ECF (km)";
    Real p_stn_top[3] "Station coordinates in topocentric (km)";
    Real N_phi "Ellipsoidal Radius of Curvature in the meridian (km)";
    Real TM[3, 3] "Transform matrix from ECF to topo";
    Integer seq[3] "rotation matrix sequence";
    Real ang[3] "sequence of rotation angles";
  equation
    N_phi = Re / sqrt(1 - ecc_e_sq * sin(d2r * stn_lat) ^ 2);
    ang = {d2r * stn_long, d2r * (90 - stn_lat), d2r * 90};
    seq = {3, 2, 3};
    p_stn_ECF[1] = (N_phi + stn_elev * 10 ^ (-3)) * cos(d2r * stn_lat) * cos(d2r * stn_long);
    p_stn_ECF[2] = (N_phi + stn_elev * 10 ^ (-3)) * cos(d2r * stn_lat) * sin(d2r * stn_long);
    p_stn_ECF[3] = (1 - ecc_e_sq) * N_phi * sin(d2r * stn_lat);
    TM = axesRotations(seq, ang);
    p_stn_top = resolve2(TM, p_stn_ECF);
  end Station;

  function Sat_ECI
    //This function will convert perifocal coordinates of Satellite to ECI
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axesRotations;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    input Real ang[3] "-argper, -inc, -raan (rad)";
    input Real p_pf[3] "Posn vector in Perifocal coords (km)";
    input Real v_pf[3] "Velocity vector in Perifocal coords (km/s)";
    output Real p_ECI[3] "Posn vector in ECI coords (km)";
    output Real v_ECI[3] "Velocity vector in ECI coords (km/s)";
  protected
    Integer seq[3] = {3, 1, 3} "Angle sequence from pf to ECI";
    Real TM[3, 3] = axesRotations(sequence = seq, angles = ang);
  algorithm
    p_ECI := resolve2(TM, p_pf);
    v_ECI := resolve2(TM, v_pf);
  end Sat_ECI;

  function theta_d
    //Calculates GMST angle
    input Real days "Number of days from J2000 to start of epoch day";
    input Real hours "day fraction from midnight to simulation time";
    output Real GMST "GMST angle (deg)";
  protected
    Real T = days / 36525 "number of julian centuries";
    Real r = 1.002737909350795 + 5.9006e-11 * T - 5.9e-15 * T ^ 2;
  algorithm
    GMST := 99.9686184 + 36000.77005 * T + 0.00038793 * T ^ 2 - 2.6e-8 * T ^ 3 "GMST at midnight [degrees]";
    GMST := GMST + 360 * r * hours;
    GMST := mod(GMST, 360) "remove integer multiples of 360 degrees";
  end theta_d;

  function sat_ECF "Converts ECI to ECF coordinates"
    // In this function, the axisRotation() function creates a 3x3 rotation matrix which is then used in the resolve2() function to rotate the p_ECI and v_ECI vectors about the inertial z-axis
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.axisRotation;
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    input Real ang "GMST angle (deg)";
    input Real p_ECI[3] "Position vector in ECI coordinates (km)";
    input Real v_ECI[3] "Velocity vector in ECI coordinates (km/s)";
    output Real p_ECF[3] "Position vector in ECF coordinates (km)";
    output Real v_ECF[3] "Relative Velocity vector in ECF coordinates (km/s)";
  protected
    Real d2r = Modelica.Constants.D2R;
    Integer ax = 3;
    Real wcross[3, 3] = skew({0., 0., 360. / 86164. * d2r});
    Real TM[3, 3] = axisRotation(axis = ax, angle = ang * d2r);
    Real v_ECF_inertial[3];
    Real v_ECI_rel[3];
  algorithm
// Determine the position vector in ECF
    p_ECF := resolve2(TM, p_ECI);
// Determine the position vector in ECF
// First, determine the ECF vecctor in the inertial frame
    v_ECF_inertial := resolve2(TM, v_ECI);
//inertial velocity vector in the ECF frame
// Second, subtract the velocity the satellite would have if it were at rest relative to the ECF at its current position
    v_ECI_rel := v_ECI - wcross * p_ECI;
    v_ECF := resolve2(TM, v_ECI_rel);
  end sat_ECF;

  function range_ECF2topo
    import Modelica.Mechanics.MultiBody.Frames.TransformationMatrices.resolve2;
    input Real p_stn_ECF[3] "Position of station in ECF coords";
    input Real p_sat_ECF[3] "Position of satellite in ECF coords";
    input Real v_sat_ECF[3] "Relative Velocity of satellite in ECF coords";
    input Real TM[3, 3] "Transform matrix from ECF to topo";
    output Real p_sat_topo[3] "Position of satellite relative to station topo coords (km)";
    output Real v_sat_topo[3] "Velocity of satellite relative to station topo coords (km/s)";
  protected
    Real p_sat_rel[3];
  algorithm
    p_sat_rel := p_sat_ECF - p_stn_ECF;
// Determine the position vector in topo
    p_sat_topo := resolve2(TM, p_sat_rel);
// Determine the velocity vector in topo
    v_sat_topo := resolve2(TM, v_sat_ECF);
  end range_ECF2topo;

  // This function calculates look angles, azimuth and elevation, at the station position in the topocentric system
  // Inputs: Topocentric position and velocity vectors of the satellite
  //Outputs: Azimuth, Elevation Look-Angles at Station

  function range_topo2look_angles
    import SI = Modelica.SIunits;
    import Modelica.SIunits.Conversions.*;
    input Real p_sat_topo[3] "Position of satellite in topo coords (km)";
    input Real v_sat_topo[3] "Velocity of satellite in topo coords (km)";
    output Real Azimuth "Azimuth look angle (deg)";
    output Real Elevation "Elevation look angle (deg)";
    output Real Azdot "Azimuth rate of change (deg/s)";
    output Real Eldot "Elevation rate of change (deg/s)";
  protected
    Real Rxy;
    Real R;
    Real pi = Modelica.Constants.pi;
  algorithm
    Azimuth := to_deg(atan2(p_sat_topo[1], p_sat_topo[2]));
    Elevation := to_deg(atan2(p_sat_topo[3], sqrt(p_sat_topo[1] ^ 2 + p_sat_topo[2] ^ 2)));
    Azimuth := mod(Azimuth, 360);
    Rxy := sqrt(p_sat_topo[1] ^ 2 + p_sat_topo[2] ^ 2);
    R := sqrt(p_sat_topo[1] ^ 2 + p_sat_topo[2] ^ 2 + p_sat_topo[3] ^ 2);
    Azdot := v_sat_topo[1] * 180 / (Rxy * pi);
    Eldot := v_sat_topo[3] * 180 * 60 / (R * pi);
  end range_topo2look_angles;
end SatTrak;