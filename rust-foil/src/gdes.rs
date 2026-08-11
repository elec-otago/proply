//! Geometry design utilities (port of `m_xgdes.f90`).
//!
//! Only `abcopy` is needed by the analysis path.

use crate::geom::lefind;
use crate::panel::{apcalc, ncalc};
use crate::spline::{scalc, segspl, seval};
use crate::s_xfoil::tecalc;
use crate::state::{IQX, Xfoil};

/// Copies the buffer airfoil to the current airfoil.
pub fn abcopy(xf: &mut Xfoil, lconf: bool) {
    if xf.nb <= 1 {
        if xf.show_output {
            eprintln!("ABCOPY: Buffer airfoil not available.");
        }
        return;
    } else if xf.nb > IQX - 5 {
        if xf.show_output {
            eprintln!("Maximum number of panel nodes  : {}", IQX - 5);
            eprintln!("Number of buffer airfoil points: {}", xf.nb);
            eprintln!("Current airfoil cannot be set.");
            eprintln!("Try executing PANE at Top Level instead.");
        }
        return;
    }

    if xf.n != xf.nb {
        xf.lbli_ni = false;
    }

    xf.n = xf.nb;
    for i in 0..xf.n {
        xf.x[i] = xf.xb[i];
        xf.y[i] = xf.yb[i];
    }
    xf.lgsame = true;

    if xf.lbflap {
        xf.xof = xf.xbf;
        xf.yof = xf.ybf;
        xf.lflap = true;
    }

    // strip out doubled points
    let mut i = 1usize;
    loop {
        if xf.x[i - 1] == xf.x[i] && xf.y[i - 1] == xf.y[i] {
            for j in i..xf.n - 1 {
                xf.x[j] = xf.x[j + 1];
                xf.y[j] = xf.y[j + 1];
            }
            xf.n -= 1;
        }
        if i >= xf.n {
            scalc(&xf.x[..xf.n], &xf.y[..xf.n], &mut xf.s[..xf.n]);
            segspl(&xf.x[..xf.n], &mut xf.xp[..xf.n], &xf.s[..xf.n]);
            segspl(&xf.y[..xf.n], &mut xf.yp[..xf.n], &xf.s[..xf.n]);

            ncalc(&xf.x[..xf.n], &xf.y[..xf.n], &xf.s[..xf.n], &mut xf.nx[..xf.n], &mut xf.ny[..xf.n]);

            lefind(
                &mut xf.sle,
                &xf.x[..xf.n],
                &xf.xp[..xf.n],
                &xf.y[..xf.n],
                &xf.yp[..xf.n],
                &xf.s[..xf.n],
                xf.show_output,
            );
            xf.xle = seval(xf.sle, &xf.x[..xf.n], &xf.xp[..xf.n], &xf.s[..xf.n]);
            xf.yle = seval(xf.sle, &xf.y[..xf.n], &xf.yp[..xf.n], &xf.s[..xf.n]);
            xf.xte = 0.5 * (xf.x[0] + xf.x[xf.n - 1]);
            xf.yte = 0.5 * (xf.y[0] + xf.y[xf.n - 1]);
            xf.chord = ((xf.xte - xf.xle).powi(2) + (xf.yte - xf.yle).powi(2)).sqrt();

            tecalc(xf);
            apcalc(xf);

            xf.lgamu = false;
            xf.lqinu = false;
            xf.lwake = false;
            xf.lqaij = false;
            xf.ladij = false;
            xf.lwdij = false;
            xf.lipan = false;
            xf.lvconv = false;
            xf.lscini = false;

            if lconf && xf.show_output {
                eprintln!(" Current airfoil nodes set from buffer airfoil nodes ( {:4} )", xf.n);
            }
            break;
        }
        i += 1;
    }
}
