// Copyright (c) Tim Molteno tim@elec.ac.nz 2026
//! Global state for the XFOIL port.
//!
//! This is the Rust translation of the Fortran module `i_xfoil` (plus the
//! parameter modules `i_pindex`, `i_blpar`, `i_circle`, `i_xbl`) and the
//! default-value initialization performed by `init` in `m_xfoil.f90`.

// ---- Primary dimensioning limit parameters ---------------------------------
pub const IQX: usize = 640; // number of surface panel nodes + 6
pub const IPX: usize = 5; // number of Qspec(s) distributions
pub const ISX: usize = 2; // number of airfoil sides
pub const IWX: usize = IQX / 8 + 2; // number of wake panel nodes
pub const IBX: usize = 4 * IQX; // number of buffer airfoil nodes
pub const IZX: usize = IQX + IWX; // number of panel nodes (airfoil + wake)
pub const IVX: usize = IQX / 2 + IWX + 50; // number of nodes along BL on one side
pub const NAX: usize = 800; // number of points in stored polar
pub const NPX: usize = 12; // number of polars and reference polars
pub const NFX: usize = 128; // number of points in one reference polar
pub const NTX: usize = 2 * IBX; // number of points in thickness/camber arrays

pub const ICX: usize = 257; // circle plane array size
pub const IMX: usize = (ICX - 1) / 4; // number of complex mapping coefficients

// ---- polar variable indexing parameters (i_pindex.f90) ---------------------
pub const IAL: usize = 1; // alpha
pub const ICL: usize = 2; // CL
pub const ICD: usize = 3; // CD
pub const ICM: usize = 4; // CM
pub const ICW: usize = 5; // CDwave
pub const ICV: usize = 6; // CDvisc
pub const ICP: usize = 7; // CDpres
pub const IMA: usize = 8; // Mach
pub const IRE: usize = 9; // Re
pub const ICH: usize = 10; // Chinge
pub const IMC: usize = 11; // Cpmin
pub const ICDH: usize = 12; // CDh
pub const ICMDOT: usize = 13; // Cmdot
pub const IPTOT: usize = 13;
pub const JNC: usize = 1; // Ncrit
pub const JTP: usize = 2; // Xtrip
pub const JTN: usize = 3; // Xtr
pub const JTI: usize = 4; // Itr
pub const JPTOT: usize = 4;

pub const NCOM: usize = 73; // number of BL Newton system unknowns per station

// ---- fundamental constants (set by init in m_xfoil.f90) ------------------
pub const PI: f64 = std::f64::consts::PI;
pub const HOPI: f64 = 0.5 / PI;
pub const QOPI: f64 = 0.25 / PI;
pub const DTOR: f64 = PI / 180.0;

// ---- state ----------------------------------------------------------------
#[derive(Clone)]
pub struct Xfoil {
    // ----- control/output switches -----
    pub show_output: bool,

    // ----- freestream / operating conditions -----
    pub adeg: f64, // angle of attack in degrees
    pub alfa: f64, // angle of attack in radians
    pub minf: f64, // freestream Mach number at current CL
    pub minf1: f64, // freestream Mach number at CL = 1
    pub minf_cl: f64, // dMINF/dCL
    pub qinf: f64, // freestream speed (defined as 1)
    pub reinf: f64, // Reynolds number for current CL
    pub reinf1: f64, // Reynolds number at CL = 1
    pub reinf_cl: f64, // dREINF/dCL
    pub clspec: f64, // specified CL
    pub retyp: i32, // Re variation with CL type
    pub matyp: i32, // Ma variation with CL type
    pub acrit: [f64; ISX], // log(critical amplification ratio)
    pub xstrip: [f64; ISX], // transition trip x/c locations
    pub tindex: [f64; ISX],
    pub waklen: f64, // wake length to chord ratio

    // ----- gas constants -----
    pub gamma: f64, // Cp/Cv
    pub gamm1: f64, // Cp/Cv - 1

    // ----- current operating point results -----
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    pub cdf: f64, // friction CD
    pub cdp: f64, // pressure CD
    pub cl_alf: f64, // dCL/dALFA
    pub cl_msq: f64, // dCL/d(MINF^2)
    pub cpmin: f64, // min Cp
    pub cpmn: f64, // min Cp (from fcpmin)
    pub cpmni: f64, // min Cp inviscid
    pub cpmnv: f64, // min Cp viscous
    pub xcpmni: f64,
    pub xcpmnv: f64,
    pub cpstar: f64, // sonic pressure coefficient
    pub qstar: f64, // sonic speed
    pub tklam: f64, // Karman-Tsien parameter
    pub tkl_msq: f64, // d(TKLAM)/d(MINF^2)

    // ----- alpha/CL setting state -----
    pub lalfa: bool, // true if alpha is specified, false if CL is specified
    pub lvisc: bool, // true if viscous option is invoked
    pub lwake: bool, // true if wake geometry has been calculated
    pub lvconv: bool, // true if converged BL solution exists
    pub lbli_ni: bool, // true if BL has been initialized (LBLini)
    pub abort_on_nan: bool, // stop iterations if NaN appears

    // ----- current airfoil geometry (panel nodes) -----
    pub n: usize, // number of points on airfoil
    pub nw: usize, // number of points in wake
    pub x: Vec<f64>, // IZX
    pub y: Vec<f64>, // IZX
    pub xp: Vec<f64>, // dX/dS
    pub yp: Vec<f64>, // dY/dS
    pub s: Vec<f64>, // arc length along airfoil
    pub sle: f64, // value of S at leading edge
    pub xle: f64, // leading edge coordinates
    pub yle: f64,
    pub xte: f64, // trailing edge coordinates
    pub yte: f64,
    pub chord: f64, // chord length
    pub cosa: f64, // cos(ALFA)
    pub sina: f64, // sin(ALFA)
    pub nx: Vec<f64>, // normal unit vector components
    pub ny: Vec<f64>,
    pub apanel: Vec<f64>, // panel angle array
    // Thickness of the "dead air" region in the wake, indexed 0-based as
    // wgap[0..nw-1].  The first wake BL station (ibl == iblte[is] + 1) reads
    // wgap[0], the last (ibl == iblte[is] + nw) reads wgap[nw-1]; in general
    // station `iw` (1-based, iw = ibl - iblte[is]) reads wgap[iw-1].  Sized to
    // IWX so that any valid iw in 1..=nw (= 1..=IWX) stays in bounds.
    pub wgap: Vec<f64>, // thickness of "dead air" region in wake (length IWX, 0-based)

    // ----- panel strengths / solution -----
    pub gam: Vec<f64>, // surface vortex panel strength array (IQX)
    pub gam_a: Vec<f64>, // dGAM/dALFA
    pub gamu: [Vec<f64>; 2], // GAM for alpha = 0, 90 deg
    pub sig: Vec<f64>, // surface/wake mass defect array (IZX)
    pub sigte: f64, // source panel strength across finite-thickness TE
    pub sigte_a: f64,
    pub gamte: f64, // vortex panel strength across finite-thickness TE
    pub gamte_a: f64,
    pub dste: f64, // TE panel length
    pub ante: f64, // projected TE thickness perp. to TE bisector
    pub aste: f64, // projected TE thickness para. to TE bisector
    pub sharp: bool, // true if DSTE == 0
    pub circ: f64, // circulation
    pub psio: f64, // streamfunction inside airfoil
    pub sst: f64, // S value at stagnation point
    pub sst_go: f64, // dSST/dGAM(IST)
    pub sst_gp: f64, // dSST/dGAM(IST+1)
    pub ist: usize, // stagnation point lies between S(IST), S(IST+1)

    // ----- inviscid influence coefficient arrays -----
    pub qinv: Vec<f64>, // tangential velocity due to surface vorticity (IZX)
    pub qinv_a: Vec<f64>, // dQINV/dalpha
    pub qvis: Vec<f64>, // tangential velocity due to vorticity & sources (IZX)
    pub qinvu: [Vec<f64>; 2], // QINV for alpha = 0, 90 deg (IZX)
    pub cpi: Vec<f64>, // inviscid Cp (IZX)
    pub cpv: Vec<f64>, // viscous Cp (IZX)
    pub qtan1: f64, // Qtan at alpha = 0 deg
    pub qtan2: f64, // Qtan at alpha = 90 deg
    pub aij: Vec<f64>, // dPsi/dGam influence matrix (IQX x IQX)
    pub aijpiv: Vec<i32>, // pivot index array
    pub bij: Vec<f64>, // dGam/dSig influence matrix (IQX x IZX)
    pub cij: Vec<f64>, // dQtan/dGam influence matrix (IWX x IQX)
    pub dij: Vec<f64>, // dQtan/dSig influence matrix (IZX x IZX)
    pub dij_t: Vec<f64>, // transposed copy of dij (row-major over j), for contiguous BL reads
    pub q: Vec<f64>, // generic coefficient matrix (IQX x IQX)
    pub dq: Vec<f64>, // generic matrix righthand side (IQX)
    pub dqdg: Vec<f64>, // dQtan/dGam (IQX)
    pub dqdm: Vec<f64>, // dQtan/dSig (IZX)
    pub dzdg: Vec<f64>, // dPsi/dGam (IQX)
    pub dzdn: Vec<f64>, // dPsi/dn (IZX)
    pub dzdm: Vec<f64>, // dPsi/dSig (IZX)
    pub z_qinf: f64, // dPsi/dQinf
    pub z_alfa: f64, // dPsi/dalfa
    pub z_qdof0: f64,
    pub z_qdof1: f64,
    pub z_qdof2: f64,
    pub z_qdof3: f64,
    pub qf0: Vec<f64>, // shape functions for QSPEC modification (IQX)
    pub qf1: Vec<f64>,
    pub qf2: Vec<f64>,
    pub qf3: Vec<f64>,
    pub qdof0: f64,
    pub qdof1: f64,
    pub qdof2: f64,
    pub qdof3: f64,
    pub lgamu: bool, // GAMU arrays exist
    pub lqinu: bool, // QINVU arrays exist
    pub lqaij: bool, // dPsi/dGam matrix computed and factored
    pub ladij: bool, // dQ/dSig matrix for airfoil computed
    pub lwdij: bool, // dQ/dSig matrix for wake computed

    // ----- matrix solve flags / params -----
    pub lipan: bool, // BL->panel pointers IPAN have been calculated
    pub itmax: i32, // max number of Newton iterations
    pub nseqex: i32, // max number of unconverged sequence points for early exit
    pub rlx: f64, // underrelaxation factor for Newton update
    pub vaccel: f64, // BL Newton acceleration parameter
    pub idamp: i32, // e^n envelope model flag

    // ----- wake -----
    pub awake: f64, // angle of attack corresponding to wake geometry (rad)
    pub avisc: f64, // angle of attack corresponding to BL solution (rad)
    pub mvisc: f64, // Mach number corresponding to BL solution
    pub algam: f64, // alpha corresponding to QGAMM distribution
    pub clgam: f64,
    pub cmgam: f64,
    pub qgamm: Vec<f64>, // surface velocity for current airfoil geometry (IBX)

    // ----- BL arrays (IVX x ISX), indexed [side][iv] -----
    pub xssi: [Vec<f64>; ISX], // BL arc length coordinate array
    pub uedg: [Vec<f64>; ISX], // BL edge velocity array
    pub uinv: [Vec<f64>; ISX], // BL edge velocity without mass defect influence
    pub uinv_a: [Vec<f64>; ISX], // dUINV/dalfa
    pub mass: [Vec<f64>; ISX], // mass defect array (= UEDG*DSTR)
    pub thet: [Vec<f64>; ISX], // momentum thickness array
    pub dstr: [Vec<f64>; ISX], // displacement thickness array
    pub tstr: [Vec<f64>; ISX], // kin. energy thickness array
    pub ctau: [Vec<f64>; ISX], // sqrt(max shear coeff) / log(amplification)
    pub tau: [Vec<f64>; ISX], // wall shear stress (plotting)
    pub dis: [Vec<f64>; ISX], // dissipation (plotting)
    pub delt: [Vec<f64>; ISX], // BL thickness (plotting)
    pub ctq: [Vec<f64>; ISX], // sqrt(equilibrium max shear coeff)
    pub vti: [Vec<f64>; ISX], // +/-1 conversion factor between panel and BL vars
    pub uslp: [Vec<f64>; ISX],
    pub guxd: [Vec<f64>; ISX], // dUe/dxBL
    pub guxq: [Vec<f64>; ISX],
    pub iblte: [i32; ISX], // BL array index at trailing edge
    pub nbl: [i32; ISX], // max BL array index
    pub ipan: [Vec<i32>; ISX], // panel index corresponding to BL location
    pub isys: [Vec<i32>; ISX], // BL Newton system line number corresponding to BL location
    pub itran: [i32; ISX], // BL array index at transition
    pub tforce: [bool; ISX], // true if transition is forced due to trip
    pub xoctr: [f64; ISX], // actual transition x/c locations
    pub yoctr: [f64; ISX], // actual transition y/c locations
    pub xssitr: [f64; ISX], // actual transition xi locations

    // ----- BL secondary variables at the "1" and "2" stations (i_xbl) -----
    pub com1: BlVars,
    pub com2: BlVars,
    pub c1sav: BlVars,
    pub c2sav: BlVars,
    pub amcrit: f64,
    pub bule: f64,
    pub cfm: f64,
    pub cfm_d1: f64,
    pub cfm_d2: f64,
    pub cfm_ms: f64,
    pub cfm_re: f64,
    pub cfm_t1: f64,
    pub cfm_t2: f64,
    pub cfm_u1: f64,
    pub cfm_u2: f64,
    pub dwte: f64,
    pub gambl: f64,
    pub gm1bl: f64,
    pub hstinv: f64,
    pub hstinv_ms: f64,
    pub hvrat: f64,
    pub qinfbl: f64,
    pub reybl: f64,
    pub reybl_ms: f64,
    pub reybl_re: f64,
    pub rstbl: f64,
    pub rstbl_ms: f64,
    pub tkbl: f64,
    pub tkbl_ms: f64,
    pub xiforc: f64,
    pub xt: f64,
    pub xt_a1: f64,
    pub xt_d1: f64,
    pub xt_d2: f64,
    pub xt_ms: f64,
    pub xt_re: f64,
    pub xt_t1: f64,
    pub xt_t2: f64,
    pub xt_u1: f64,
    pub xt_u2: f64,
    pub xt_x1: f64,
    pub xt_x2: f64,
    pub xt_xf: f64,
    pub idampv: i32,
    pub simi: bool,
    pub tran: bool,
    pub trforc: bool,
    pub trfree: bool,
    pub turb: bool,
    pub wake: bool,
    pub vs1: [[f64; 5]; 4],
    pub vs2: [[f64; 5]; 4],
    pub vsm: [f64; 4],
    pub vsr: [f64; 4],
    pub vsrez: [f64; 4],
    pub vsx: [f64; 4],

    // ----- BL Newton system arrays -----
    pub nsys: usize, // total number of lines in BL Newton system
    pub va: Vec<f64>, // diagonal blocks (3 x 2 x IZX)
    pub vb: Vec<f64>, // off-diagonal blocks (3 x 2 x IZX)
    pub vz: Vec<f64>, // way-off-diagonal block at TE station (3 x 2)
    pub vm: Vec<f64>, // mass-influence coefficient vectors (3 x IZX x IZX)
    pub vdel: Vec<f64>, // residual and solution vectors (3 x 2 x IZX)
    pub rmsbl: f64, // rms change from BL Newton system solution
    pub rmxbl: f64, // max change from BL Newton system solution
    pub imxbl: i32, // location of max change
    pub ismxbl: i32, // index of BL side containing max change
    pub vmxbl: char, // character identifying variable with max change

    // ----- work arrays -----
    pub w1: Vec<f64>, // 6*IQX
    pub w2: Vec<f64>,
    pub w3: Vec<f64>,
    pub w4: Vec<f64>,
    pub w5: Vec<f64>,
    pub w6: Vec<f64>,
    pub w7: Vec<f64>,
    pub w8: Vec<f64>,

    // ----- BL Newton scratch arrays (reused across iterations; see setbl/update) -----
    pub bl_usav: [Vec<f64>; 2],   // [2][IVX]
    pub bl_ule1_m: Vec<f64>,      // 2*IVX
    pub bl_ule2_m: Vec<f64>,      // 2*IVX
    pub bl_ute1_m: Vec<f64>,      // 2*IVX
    pub bl_ute2_m: Vec<f64>,      // 2*IVX
    pub bl_u1_m: Vec<f64>,        // 2*IVX
    pub bl_d1_m: Vec<f64>,        // 2*IVX
    pub bl_u2_m: Vec<f64>,        // 2*IVX
    pub bl_d2_m: Vec<f64>,        // 2*IVX
    pub bl_unew: [Vec<f64>; 2],   // [2][IVX]
    pub bl_u_ac: [Vec<f64>; 2],   // [2][IVX]
    pub bl_qnew: Vec<f64>,        // IQX
    pub bl_q_ac: Vec<f64>,        // IQX

    // ----- buffer airfoil -----
    pub nb: usize, // number of points in buffer airfoil array
    pub xb: Vec<f64>, // buffer airfoil coordinate arrays (IBX)
    pub yb: Vec<f64>,
    pub xbp: Vec<f64>, // dXB/dSB
    pub ybp: Vec<f64>, // dYB/dSB
    pub sb: Vec<f64>, // spline parameter for buffer airfoil
    pub snew: Vec<f64>, // new panel endpoint arc length array (5*IBX)
    pub sble: f64, // LE tangency-point SB location
    pub chordb: f64, // chord
    pub areab: f64, // area
    pub radble: f64, // LE radius
    pub angbte: f64, // TE angle (rad)
    pub ei11ba: f64,
    pub ei22ba: f64,
    pub apx1ba: f64,
    pub apx2ba: f64,
    pub ei11bt: f64,
    pub ei22bt: f64,
    pub apx1bt: f64,
    pub apx2bt: f64,
    pub thickb: f64, // max thickness
    pub cambrb: f64, // max camber
    pub xbf: f64, // buffer airfoil flap hinge coordinates
    pub ybf: f64,
    pub lbflap: bool,
    pub lnorm: bool, // normalize input buffer airfoil
    pub lclock: bool, // source airfoil coordinates are clockwise
    pub lgsame: bool, // current and buffer airfoils are identical

    // ----- current airfoil flap -----
    pub xof: f64,
    pub yof: f64,
    pub hmom: f64,
    pub hfx: f64,
    pub hfy: f64,
    pub lflap: bool,

    // ----- paneling parameters -----
    pub npan: usize, // default/specified number of points on airfoil
    pub cvpar: f64, // curvature attraction parameter
    pub cterat: f64, // TE panel density / LE panel density ratio
    pub ctrrat: f64, // local refinement panel density / LE panel density ratio
    pub xsref1: f64, // suction side local refinement x/c limits
    pub xsref2: f64,
    pub xpref1: f64, // pressure side local refinement x/c limits
    pub xpref2: f64,

    // ----- moment reference -----
    pub xcmref: f64,
    pub ycmref: f64,

    // ----- circle plane state (i_circle.f90) -----
    pub nc1: usize, // number of circle plane points
    pub mc: i32,
    pub mct: i32,
    pub nc: i32,
    pub ag0: f64,
    pub agte: f64,
    pub dwc: f64,
    pub qim0: f64,
    pub qimold: f64,
    pub zleold: num_complex::Complex64,
    pub chordz: num_complex::Complex64,
    pub dzte: num_complex::Complex64,
    pub cn: Vec<num_complex::Complex64>, // IMX+1
    pub piq: Vec<num_complex::Complex64>, // ICX
    pub zc: Vec<num_complex::Complex64>, // ICX
    pub zcoldw: Vec<num_complex::Complex64>, // ICX
    pub zc_cn: Vec<num_complex::Complex64>, // ICX * IMX/4
    pub eiw: Vec<num_complex::Complex64>, // ICX * (IMX+1)
    pub sc: Vec<f64>, // ICX
    pub scold: Vec<f64>,
    pub wc: Vec<f64>,
    pub xcold: Vec<f64>,
    pub ycold: Vec<f64>,
    pub leiw: bool, // unit circle complex array initialized
    pub lscini: bool, // old-airfoil circle-plane arc length s(w) exists

    // ----- BL parameter calibration constants (i_blpar.f90) -----
    pub cffac: f64,
    pub ctcon: f64,
    pub ctrcex: f64,
    pub ctrcon: f64,
    pub dlcon: f64,
    pub duxcon: f64,
    pub gacon: f64,
    pub gbcon: f64,
    pub gccon: f64,
    pub sccon: f64,

    // ----- misc -----
    pub kimage: i32,
    pub yimage: f64,
    pub limage: bool,
    pub kdelim: i32,
    pub xcmax: f64,
    pub xcmin: f64,
    pub ycmax: f64,
    pub ycmin: f64,
    pub ncpref: i32,
    pub nname: usize,
    pub name: String,
    pub nprefix: usize,
    pub prefix: String,
    pub oname: String,
    pub pfname: [String; NPX],
    pub ispars: String,
}

impl Default for Xfoil {
    fn default() -> Self {
        Self::new()
    }
}

impl Xfoil {
    pub fn new() -> Self {
        let mut xf = Xfoil {
            show_output: true,

            adeg: 0.0,
            alfa: 0.0,
            minf: 0.0,
            minf1: 0.0,
            minf_cl: 0.0,
            qinf: 1.0,
            reinf: 0.0,
            reinf1: 0.0,
            reinf_cl: 0.0,
            clspec: 0.0,
            retyp: 1,
            matyp: 1,
            acrit: [9.0; ISX],
            xstrip: [1.0; ISX],
            tindex: [0.0; ISX],
            waklen: 1.0,

            gamma: 1.4,
            gamm1: 0.4,

            cl: 0.0,
            cm: 0.0,
            cd: 0.0,
            cdf: 0.0,
            cdp: 0.0,
            cl_alf: 0.0,
            cl_msq: 0.0,
            cpmin: 0.0,
            cpmn: 0.0,
            cpmni: 0.0,
            cpmnv: 0.0,
            xcpmni: 0.0,
            xcpmnv: 0.0,
            cpstar: 0.0,
            qstar: 0.0,
            tklam: 0.0,
            tkl_msq: 0.0,

            lalfa: true,
            lvisc: false,
            lwake: false,
            lvconv: false,
            lbli_ni: false,
            abort_on_nan: true,

            n: 0,
            nw: 0,
            x: vec![0.0; IZX],
            y: vec![0.0; IZX],
            xp: vec![0.0; IZX],
            yp: vec![0.0; IZX],
            s: vec![0.0; IZX],
            sle: 0.0,
            xle: 0.0,
            yle: 0.0,
            xte: 0.0,
            yte: 0.0,
            chord: 0.0,
            cosa: 1.0,
            sina: 0.0,
            nx: vec![0.0; IZX],
            ny: vec![0.0; IZX],
            apanel: vec![0.0; IZX],
            wgap: vec![0.0; IWX],

            gam: vec![0.0; IQX],
            gam_a: vec![0.0; IQX],
            gamu: [vec![0.0; IQX], vec![0.0; IQX]],
            sig: vec![0.0; IZX],
            sigte: 0.0,
            sigte_a: 0.0,
            gamte: 0.0,
            gamte_a: 0.0,
            dste: 0.0,
            ante: 0.0,
            aste: 0.0,
            sharp: false,
            circ: 0.0,
            psio: 0.0,
            sst: 0.0,
            sst_go: 0.0,
            sst_gp: 0.0,
            ist: 0,

            qinv: vec![0.0; IZX],
            qinv_a: vec![0.0; IZX],
            qvis: vec![0.0; IZX],
            qinvu: [vec![0.0; IZX], vec![0.0; IZX]],
            cpi: vec![0.0; IZX],
            cpv: vec![0.0; IZX],
            qtan1: 0.0,
            qtan2: 0.0,
            aij: vec![0.0; IQX * IQX],
            aijpiv: vec![0; IQX],
            bij: vec![0.0; IQX * IZX],
            cij: vec![0.0; IWX * IQX],
            dij: vec![0.0; IZX * IZX],
            dij_t: vec![0.0; IZX * IZX],
            q: vec![0.0; IQX * IQX],
            dq: vec![0.0; IQX],
            dqdg: vec![0.0; IQX],
            dqdm: vec![0.0; IZX],
            dzdg: vec![0.0; IQX],
            dzdn: vec![0.0; IZX],
            dzdm: vec![0.0; IZX],
            z_qinf: 0.0,
            z_alfa: 0.0,
            z_qdof0: 0.0,
            z_qdof1: 0.0,
            z_qdof2: 0.0,
            z_qdof3: 0.0,
            qf0: vec![0.0; IQX],
            qf1: vec![0.0; IQX],
            qf2: vec![0.0; IQX],
            qf3: vec![0.0; IQX],
            qdof0: 0.0,
            qdof1: 0.0,
            qdof2: 0.0,
            qdof3: 0.0,
            lgamu: false,
            lqinu: false,
            lqaij: false,
            ladij: false,
            lwdij: false,

            lipan: false,
            itmax: 20,
            nseqex: 4,
            rlx: 0.0,
            vaccel: 0.01,
            idamp: 0,

            awake: 0.0,
            avisc: 0.0,
            mvisc: 0.0,
            algam: 0.0,
            clgam: 0.0,
            cmgam: 0.0,
            qgamm: vec![0.0; IBX],

            xssi: [vec![0.0; IVX], vec![0.0; IVX]],
            uedg: [vec![0.0; IVX], vec![0.0; IVX]],
            uinv: [vec![0.0; IVX], vec![0.0; IVX]],
            uinv_a: [vec![0.0; IVX], vec![0.0; IVX]],
            mass: [vec![0.0; IVX], vec![0.0; IVX]],
            thet: [vec![0.0; IVX], vec![0.0; IVX]],
            dstr: [vec![0.0; IVX], vec![0.0; IVX]],
            tstr: [vec![0.0; IVX], vec![0.0; IVX]],
            ctau: [vec![0.0; IVX], vec![0.0; IVX]],
            tau: [vec![0.0; IVX], vec![0.0; IVX]],
            dis: [vec![0.0; IVX], vec![0.0; IVX]],
            delt: [vec![0.0; IVX], vec![0.0; IVX]],
            ctq: [vec![0.0; IVX], vec![0.0; IVX]],
            vti: [vec![0.0; IVX], vec![0.0; IVX]],
            uslp: [vec![0.0; IVX], vec![0.0; IVX]],
            guxd: [vec![0.0; IVX], vec![0.0; IVX]],
            guxq: [vec![0.0; IVX], vec![0.0; IVX]],
            iblte: [0; ISX],
            nbl: [0; ISX],
            ipan: [vec![0; IVX], vec![0; IVX]],
            isys: [vec![0; IVX], vec![0; IVX]],
            itran: [0; ISX],
            tforce: [false; ISX],
            xoctr: [1.0; ISX],
            yoctr: [0.0; ISX],
            xssitr: [0.0; ISX],

            com1: BlVars::default(),
            com2: BlVars::default(),
            c1sav: BlVars::default(),
            c2sav: BlVars::default(),
            amcrit: 0.0,
            bule: 0.0,
            cfm: 0.0,
            cfm_d1: 0.0,
            cfm_d2: 0.0,
            cfm_ms: 0.0,
            cfm_re: 0.0,
            cfm_t1: 0.0,
            cfm_t2: 0.0,
            cfm_u1: 0.0,
            cfm_u2: 0.0,
            dwte: 0.0,
            gambl: 0.0,
            gm1bl: 0.0,
            hstinv: 0.0,
            hstinv_ms: 0.0,
            // Upstream XFOIL value (HVRAIN in i_blpar.f90).  The xfoil-python
            // port dropped this DATA statement and effectively used 0.0; we use
            // the canonical XFOIL value so results match Drela's code rather
            // than the Python port.  See README "Fixes vs upstream".
            hvrat: 0.25,
            qinfbl: 0.0,
            reybl: 0.0,
            reybl_ms: 0.0,
            reybl_re: 0.0,
            rstbl: 0.0,
            rstbl_ms: 0.0,
            tkbl: 0.0,
            tkbl_ms: 0.0,
            xiforc: 0.0,
            xt: 0.0,
            xt_a1: 0.0,
            xt_d1: 0.0,
            xt_d2: 0.0,
            xt_ms: 0.0,
            xt_re: 0.0,
            xt_t1: 0.0,
            xt_t2: 0.0,
            xt_u1: 0.0,
            xt_u2: 0.0,
            xt_x1: 0.0,
            xt_x2: 0.0,
            xt_xf: 0.0,
            idampv: 0,
            simi: false,
            tran: false,
            trforc: false,
            trfree: false,
            turb: false,
            wake: false,
            vs1: [[0.0; 5]; 4],
            vs2: [[0.0; 5]; 4],
            vsm: [0.0; 4],
            vsr: [0.0; 4],
            vsrez: [0.0; 4],
            vsx: [0.0; 4],

            nsys: 0,
            va: vec![0.0; 3 * 2 * IZX],
            vb: vec![0.0; 3 * 2 * IZX],
            vz: vec![0.0; 3 * 2],
            vm: vec![0.0; 3 * IZX * IZX],
            vdel: vec![0.0; 3 * 2 * IZX],
            rmsbl: 0.0,
            rmxbl: 0.0,
            imxbl: 0,
            ismxbl: 0,
            vmxbl: ' ',

            w1: vec![0.0; 6 * IQX],
            w2: vec![0.0; 6 * IQX],
            w3: vec![0.0; 6 * IQX],
            w4: vec![0.0; 6 * IQX],
            w5: vec![0.0; 6 * IQX],
            w6: vec![0.0; 6 * IQX],
            w7: vec![0.0; 6 * IQX],
            w8: vec![0.0; 6 * IQX],

            bl_usav: [vec![0.0; IVX], vec![0.0; IVX]],
            bl_ule1_m: vec![0.0; 2 * IVX],
            bl_ule2_m: vec![0.0; 2 * IVX],
            bl_ute1_m: vec![0.0; 2 * IVX],
            bl_ute2_m: vec![0.0; 2 * IVX],
            bl_u1_m: vec![0.0; 2 * IVX],
            bl_d1_m: vec![0.0; 2 * IVX],
            bl_u2_m: vec![0.0; 2 * IVX],
            bl_d2_m: vec![0.0; 2 * IVX],
            bl_unew: [vec![0.0; IVX], vec![0.0; IVX]],
            bl_u_ac: [vec![0.0; IVX], vec![0.0; IVX]],
            bl_qnew: vec![0.0; IQX],
            bl_q_ac: vec![0.0; IQX],

            nb: 0,
            xb: vec![0.0; IBX],
            yb: vec![0.0; IBX],
            xbp: vec![0.0; IBX],
            ybp: vec![0.0; IBX],
            sb: vec![0.0; IBX],
            snew: vec![0.0; 5 * IBX],
            sble: 0.0,
            chordb: 0.0,
            areab: 0.0,
            radble: 0.0,
            angbte: 0.0,
            ei11ba: 0.0,
            ei22ba: 0.0,
            apx1ba: 0.0,
            apx2ba: 0.0,
            ei11bt: 0.0,
            ei22bt: 0.0,
            apx1bt: 0.0,
            apx2bt: 0.0,
            thickb: 0.0,
            cambrb: 0.0,
            xbf: 0.0,
            ybf: 0.0,
            lbflap: false,
            lnorm: true,
            lclock: false,
            lgsame: false,

            xof: 0.0,
            yof: 0.0,
            hmom: 0.0,
            hfx: 0.0,
            hfy: 0.0,
            lflap: false,

            npan: 160,
            cvpar: 1.0,
            cterat: 0.15,
            ctrrat: 0.2,
            xsref1: 1.0,
            xsref2: 1.0,
            xpref1: 1.0,
            xpref2: 1.0,

            xcmref: 0.25,
            ycmref: 0.0,

            nc1: 0,
            mc: 0,
            mct: 0,
            nc: 0,
            ag0: 0.0,
            agte: 0.0,
            dwc: 0.0,
            qim0: 0.0,
            qimold: 0.0,
            zleold: num_complex::Complex64::new(0.0, 0.0),
            chordz: num_complex::Complex64::new(0.0, 0.0),
            dzte: num_complex::Complex64::new(0.0, 0.0),
            cn: vec![num_complex::Complex64::new(0.0, 0.0); IMX + 1],
            piq: vec![num_complex::Complex64::new(0.0, 0.0); ICX],
            zc: vec![num_complex::Complex64::new(0.0, 0.0); ICX],
            zcoldw: vec![num_complex::Complex64::new(0.0, 0.0); ICX],
            zc_cn: vec![num_complex::Complex64::new(0.0, 0.0); ICX * (IMX / 4)],
            eiw: vec![num_complex::Complex64::new(0.0, 0.0); ICX * (IMX + 1)],
            sc: vec![0.0; ICX],
            scold: vec![0.0; ICX],
            wc: vec![0.0; ICX],
            xcold: vec![0.0; ICX],
            ycold: vec![0.0; ICX],
            leiw: false,
            lscini: false,

            cffac: 0.0,
            ctcon: 0.0,
            ctrcex: 0.0,
            ctrcon: 0.0,
            dlcon: 0.0,
            duxcon: 0.0,
            gacon: 0.0,
            gbcon: 0.0,
            gccon: 0.0,
            sccon: 0.0,

            kimage: 1,
            yimage: -10.0,
            limage: false,
            kdelim: 0,
            xcmax: 0.0,
            xcmin: 0.0,
            ycmax: 0.0,
            ycmin: 0.0,
            ncpref: 0,
            nname: 32,
            name: "                                ".to_string(),
            nprefix: 0,
            prefix: " ".to_string(),
            oname: " ".to_string(),
            pfname: std::array::from_fn(|_| " ".to_string()),
            ispars: " -2.0  3.0  -2.5  3.5".to_string(),
        };

        // ---- set unity freestream speed ----
        xf.qinf = 1.0;

        // ---- initialize freestream Mach number to zero ----
        xf.matyp = 1;
        xf.minf1 = 0.0;

        // ---- circle plane array size (largest 2^n + 1 that fits) ----
        let ann = ((2 * IQX - 1) as f64).ln() / 2.0f64.ln();
        let nn = (ann + 0.00001) as i32;
        xf.nc1 = ((1usize << nn) + 1).min(257);
        let _ = &mut xf.mc;

        // ---- set MINF, REINF, based on current CL-dependence ----
        xf.minf_cl = 0.0;
        xf.reinf_cl = 0.0;

        xf
    }

    /// Index into the flat 3x2xIZX BL block arrays.  All indices are 1-based,
    /// matching the Fortran `va(3, 2, IZX)` / `vdel(3, 2, IZX)` layout.
    #[inline]
    pub fn v_index(iv: usize, l: usize, k: usize) -> usize {
        ((iv - 1) * 2 + (l - 1)) * 3 + (k - 1)
    }

    /// Index into the flat 3xIZXxIZX mass-influence array, matching the
    /// Fortran `vm(3, IZX, IZX)` layout: element (k, j, i) with 1-based
    /// indices lives at `((i-1)*IZX + (j-1))*3 + (k-1)`.
    #[inline]
    pub fn vm_index(i: usize, j: usize, k: usize) -> usize {
        ((i - 1) * IZX + (j - 1)) * 3 + (k - 1)
    }

    /// Index into the flat IQX x IQX matrix, Fortran column-major layout:
    /// element (i, j) [0-based] of `aij(IQX, IQX)` lives at `j*IQX + i`.
    #[inline]
    pub fn m_index(i: usize, j: usize) -> usize {
        j * IQX + i
    }

    /// Index into the flat IQX x IZX matrix, Fortran column-major layout:
    /// element (i, j) [0-based] of `bij(IQX, IZX)` lives at `j*IQX + i`.
    #[inline]
    pub fn b_index(i: usize, j: usize) -> usize {
        j * IQX + i
    }

    /// Index into the flat IZX x IZX matrix, Fortran column-major layout:
    /// element (i, j) [0-based] of `dij(IZX, IZX)` lives at `j*IZX + i`.
    #[inline]
    pub fn d_index(i: usize, j: usize) -> usize {
        j * IZX + i
    }
}

/// Secondary BL variables at one station (port of the `i_xbl` named
/// pointers into `com1`/`com2`, each 73 entries long).  The field order
/// matches the pointer layout established by `preptrs` in `m_xbl.f90`, so
/// a whole-station copy via `com1 = com2` is equivalent to the Fortran
/// `do icom = 1, NCOM; COM1(icom) = COM2(icom); enddo`.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlVars {
    pub x: f64, // 1   station arc length
    pub u: f64, // 2   edge velocity
    pub t: f64, // 3   momentum thickness
    pub d: f64, // 4   displacement thickness (minus wake gap)
    pub s: f64, // 5   sqrt(shear coeff) (laminar) / shear coeff (turb)
    pub ampl: f64, // 6   log amplification ratio
    pub u_uei: f64, // 7   dU/dUei
    pub u_ms: f64, // 8   dU/dMsq
    pub dw: f64, // 9   wake gap
    pub h: f64, // 10  shape parameter D/T
    pub h_t: f64, // 11
    pub h_d: f64, // 12
    pub m: f64, // 13  edge Mach number squared
    pub m_u: f64, // 14
    pub m_ms: f64, // 15
    pub r: f64, // 16  edge density
    pub r_u: f64, // 17
    pub r_ms: f64, // 18
    pub v: f64, // 19  molecular viscosity
    pub v_u: f64, // 20
    pub v_ms: f64, // 21
    pub v_re: f64, // 22
    pub hk: f64, // 23  kinematic shape parameter
    pub hk_u: f64, // 24
    pub hk_t: f64, // 25
    pub hk_d: f64, // 26
    pub hk_ms: f64, // 27
    pub hs: f64, // 28  kinetic energy shape parameter
    pub hs_u: f64, // 29
    pub hs_t: f64, // 30
    pub hs_d: f64, // 31
    pub hs_ms: f64, // 32
    pub hs_re: f64, // 33
    pub hc: f64, // 34  density shape parameter
    pub hc_u: f64, // 35
    pub hc_t: f64, // 36
    pub hc_d: f64, // 37
    pub hc_ms: f64, // 38
    pub rt: f64, // 39  momentum thickness Reynolds number
    pub rt_u: f64, // 40
    pub rt_t: f64, // 41
    pub rt_ms: f64, // 42
    pub rt_re: f64, // 43
    pub cf: f64, // 44  skin friction coefficient
    pub cf_u: f64, // 45
    pub cf_t: f64, // 46
    pub cf_d: f64, // 47
    pub cf_ms: f64, // 48
    pub cf_re: f64, // 49
    pub di: f64, // 50  dissipation function 2CD/H*
    pub di_u: f64, // 51
    pub di_t: f64, // 52
    pub di_d: f64, // 53
    pub di_s: f64, // 54
    pub di_ms: f64, // 55
    pub di_re: f64, // 56
    pub us: f64, // 57  normalized slip velocity
    pub us_u: f64, // 58
    pub us_t: f64, // 59
    pub us_d: f64, // 60
    pub us_ms: f64, // 61
    pub us_re: f64, // 62
    pub cq: f64, // 63  equilibrium wall shear coeff ** 1/2
    pub cq_u: f64, // 64
    pub cq_t: f64, // 65
    pub cq_d: f64, // 66
    pub cq_ms: f64, // 67
    pub cq_re: f64, // 68
    pub de: f64, // 69  BL thickness
    pub de_u: f64, // 70
    pub de_t: f64, // 71
    pub de_d: f64, // 72
    pub de_ms: f64, // 73
}

// ---- complex number support (i_circle.f90) --------------------------------
pub mod num_complex {
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Complex64 {
        pub re: f64,
        pub im: f64,
    }

    impl Complex64 {
        pub fn new(re: f64, im: f64) -> Self {
            Complex64 { re, im }
        }
        pub fn norm(&self) -> f64 {
            (self.re * self.re + self.im * self.im).sqrt()
        }
        pub fn arg(&self) -> f64 {
            self.im.atan2(self.re)
        }
    }

    impl std::ops::Add for Complex64 {
        type Output = Complex64;
        fn add(self, rhs: Complex64) -> Complex64 {
            Complex64::new(self.re + rhs.re, self.im + rhs.im)
        }
    }
    impl std::ops::Sub for Complex64 {
        type Output = Complex64;
        fn sub(self, rhs: Complex64) -> Complex64 {
            Complex64::new(self.re - rhs.re, self.im - rhs.im)
        }
    }
    impl std::ops::Mul for Complex64 {
        type Output = Complex64;
        fn mul(self, rhs: Complex64) -> Complex64 {
            Complex64::new(
                self.re * rhs.re - self.im * rhs.im,
                self.re * rhs.im + self.im * rhs.re,
            )
        }
    }
    impl std::ops::Div for Complex64 {
        type Output = Complex64;
        fn div(self, rhs: Complex64) -> Complex64 {
            let d = rhs.re * rhs.re + rhs.im * rhs.im;
            Complex64::new(
                (self.re * rhs.re + self.im * rhs.im) / d,
                (self.im * rhs.re - self.re * rhs.im) / d,
            )
        }
    }
    impl std::ops::AddAssign for Complex64 {
        fn add_assign(&mut self, rhs: Complex64) {
            self.re += rhs.re;
            self.im += rhs.im;
        }
    }
    impl std::ops::MulAssign for Complex64 {
        fn mul_assign(&mut self, rhs: Complex64) {
            let re = self.re * rhs.re - self.im * rhs.im;
            let im = self.re * rhs.im + self.im * rhs.re;
            self.re = re;
            self.im = im;
        }
    }
    impl std::ops::Neg for Complex64 {
        type Output = Complex64;
        fn neg(self) -> Complex64 {
            Complex64::new(-self.re, -self.im)
        }
    }
}
