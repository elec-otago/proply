# Super simple prop design by Tim Molteno
# Copyright (c) 2016-2017. Tim Molteno tim@molteno.net
#

import math
import textwrap
import logging

import numpy as np

from proply import foil
from proply import stl_tools
from proply import motor_model
from proply import optimize
from proply.smooth import smooth

from proply.blade_element import BladeElement

from scipy.interpolate import PchipInterpolator, interp1d

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from proply import templates

logger = logging.getLogger(__name__)


class Prop:
    """
    A prop is a collection of BladeElement objects.
    """

    def __init__(self, param, resolution):
        self.param = param
        self.radial_resolution = resolution  # How often to create a profile
        self.radial_steps = int(self.param.radius / self.radial_resolution)
        self.n_blades = 2
        self.max_depth_interpolator = None
        self.scimitar_interpolator = None

    def new_blade_element(self, foilclass, r, rpm, twist):
        y_limit = self.get_max_depth(r)
        x_limit = self.get_max_chord(r, twist)
        thickness = self.get_foil_thickness(r)

        f = foilclass(chord=x_limit, thickness=thickness / x_limit)
        f.set_trailing_edge(self.param.trailing_edge / 1000.0)

        c_max = f.get_max_chord(x_limit, y_limit, twist)
        print(("Max Chord {}".format(c_max)))
        f.modify_chord(c_max)

        be = BladeElement(
            r,
            dr=self.radial_resolution,
            foil=f,
            twist=twist,
            rpm=rpm,
            u_0=self.param.forward_airspeed,
        )
        return be

    def new_foil(self, r, rpm, twist):
        return self.new_blade_element(foil.Foil, r, rpm, twist)

    def get_air_velocity_at_prop(self, torque, rpm, rho=1.225):
        r = self.param.radius
        tau = torque
        rps = rpm / 60.0
        v = np.sqrt((r ** 2 * rho) ** (1 / 3) * (rps * tau) ** (2 / 3) / (r ** 2 * rho))
        return v / 2.0  # Half of the

    def get_air_density(self):
        return 1.225  # kg m^{-3}

    def get_max_chord(self, r, twist):
        """Allowed chord as a function of radius (m)
        Limited by mechanical strength, or weight issues

        k/r_max**2 = tip_chord => k = tip_chord * r_max**2
        c = k / r**2
        """
        k = self.param.tip_chord * self.param.radius**2
        c = k / r**2

        upper_limit = (2.0 * np.pi * r / (self.n_blades + 2.0)) / np.cos(twist)

        c = min(c, upper_limit)
        return c

    def get_scimitar_offset(self, r):
        """How much forward or aft of the centerline to place the foil"""
        if self.scimitar_interpolator is None:
            hub_r = self.param.hub_radius
            max_r = self.param.radius * 0.8

            hub_c = 0.0
            max_c = self.param.radius * (self.param.scimitar_percent / 100.0)
            end_c = 0.0
            x = np.array([0, hub_r, max_r, self.param.radius])
            y = np.array([hub_c, 1.1 * hub_c, max_c, end_c])

            self.scimitar_interpolator = PchipInterpolator(x, y)

            # import matplotlib.pyplot as plt
            # rpts = np.linspace(0, self.param.radius, 40)
            # plt.plot(rpts, self.scimitar_interpolator(rpts), label='scimitar')
            # plt.plot(x, y, 'x', label='points')
            # plt.legend()
            # plt.grid(True)
            # plt.xlabel('Angle of Attack')
            # plt.show()
        return self.scimitar_interpolator(r)

    def get_foil_thickness(self, r):
        """Allowed foil thickness as a function of radius (m)
        Limited by mechanical strength, or weight issues
        """
        thickness_root = self.param.hub_depth * 1.0
        thickness_end = self.param.hub_depth * 0.1
        # Solve s + kr^3 = end && s + kh^3 = start
        # Subtract kr^3 - k h^3 = (end - start) => k = (end - start) / (r^3 - h^3)
        # s = end - kr^3
        p = 0.3
        k = (thickness_root - thickness_end) / (
            self.param.hub_radius ** p - self.param.radius ** p
        )
        s = thickness_end - k * self.param.radius ** p
        thickness = s + k * r ** (p)

        # import matplotlib.pyplot as plt
        # rpts = np.linspace(self.param.hub_radius, self.param.radius, 40)
        # plt.plot(rpts, s + k*rpts**-p, label='thickness')
        # plt.legend()
        # plt.grid(True)
        # plt.xlabel('r')
        # plt.show()

        return thickness

    def get_max_depth(self, r):
        """Allowed depth of the prop as a function of radius (m)
        This is a property of the environment that the prop operates in.

        TODO Load this from the exclude zone of the prop description
        """
        if self.max_depth_interpolator is None:
            hub_r = self.param.hub_radius
            hub_depth = self.param.hub_depth
            max_depth = self.param.hub_depth * 3
            max_r = self.param.radius / 3.0
            end_depth = self.param.hub_depth * 2

            x = np.array(
                [
                    0,
                    hub_r / 2,
                    hub_r,
                    1.1 * hub_r,
                    1.5 * hub_r,
                    max_r,
                    0.9 * self.param.radius,
                    self.param.radius,
                ]
            )
            y = np.array(
                [
                    hub_depth,
                    hub_depth,
                    hub_depth,
                    hub_depth,
                    1.1 * hub_depth,
                    max_depth,
                    1.2 * end_depth,
                    end_depth,
                ]
            )
            self.max_depth_interpolator = interp1d(x, y, "linear")

            import matplotlib.pyplot as plt

            rpts = np.linspace(0, self.param.radius, 40)
            plt.plot(rpts, self.max_depth_interpolator(rpts), label="max depth")
            plt.plot(rpts, self.get_foil_thickness(rpts), label="foil thickness")
            plt.plot(x, y, "x", label="points")

            plt.legend()
            plt.grid(True)
            plt.xlabel("r")
            plt.show()

        return self.max_depth_interpolator(r)

    def get_forces(self, rpm):
        torque = 0.0
        thrust = 0.0
        for be in self.blade_elements:
            be.rpm = rpm
            dv_goal = be.dv
            dv, a_prime, err = be.bem(self.n_blades)

            if err < 0.01:
                dT = be.dT()
                dM = be.dM()
                thrust += dT
                torque += dM
                logger.info(
                    "r={}, theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={} ".format(
                        be.r, np.degrees(be.get_twist()), dv, a_prime, dT, dM, dT / dM
                    )
                )

            else:
                logger.warning(
                    "r={}: BEM did not converge {} {} {}".format(
                        be.r, be.dv, dv_goal, a_prime, err
                    )
                )
                if err > 0.5:
                    be.dv = 1.0
                    be.a_prime = 0.0

        return torque, thrust

    def gen_mesh(self, filename, n):
        import pygmsh as pg

        geom = pg.Geometry()

        loops = []
        for be in self.blade_elements:
            car = 0.5 / 1000
            line_l, line_u = be.get_foil_points(n, self.get_scimitar_offset(be.r))
            loop_points = np.concatenate((line_l, line_u[::-1]), axis=0)
            g_pts = []
            for p in loop_points[0:-2]:
                g_pts.append(geom.add_point(p, car))

            l_foil = geom.add_bspline(g_pts)

            loops.append(l_foil)

        # print geom.get_code()
        geom.add_ruled_surface([loops[0], loops[1]])

        # l in range(len(loops) - 1):
        # geom.add_surface(loops[l], loops[l+1])
        # poly = geom.add_polygon([
        # [0.0,   0.5, 0.0],
        # [-0.1,  0.1, 0.0],
        # [-0.5,  0.0, 0.0],
        # [-0.1, -0.1, 0.0],
        # [0.0,  -0.5, 0.0],
        # [0.1,  -0.1, 0.0],
        # [0.5,   0.0, 0.0],
        # [0.1,   0.1, 0.0]
        # ],
        # lcar=0.05
        # )
        # axis = [0, 0, 1]

        # geom.extrude(
        #'Surface{%s}' % poly,
        # translation_axis=axis,
        # rotation_axis=axis,
        # point_on_axis=[0, 0, 0],
        # angle=2.0 / 6.0 * np.pi
        # )

        points, cells = pg.generate_mesh(geom)

        import meshio

        meshio.write(filename, points, cells)

    def gen_stl(self, filename, n):
        stl = stl_tools.STL()

        scale = 1000.0  # Convert to mm.
        top_lines = []
        bottom_lines = []

        bottom_edge = []
        top_edge = []

        hub_element = self.blade_elements[0]
        x0, x1, y0, y1 = hub_element.foil.get_bounding_box(hub_element.get_twist())

        for be in self.blade_elements:
            line_l, line_u = be.get_foil_points(n, self.get_scimitar_offset(be.r))

            top_lines.append(line_u * scale)
            top_edge.append(line_u[-1, :] * scale)

            bottom_lines.append(line_l * scale)
            bottom_edge.append(line_l[-1, :] * scale)

        stl.add_line(bottom_lines[0])
        ## Do the top surface
        for tl in top_lines:
            stl.add_line(tl)

        ## Do the bottom surface
        bottom_lines.reverse()
        for bl in bottom_lines:
            stl.add_line(bl)

        ## Join the trailing edges
        stl.new_block()
        stl.add_line(bottom_edge)
        stl.add_line(top_edge)

        stl.gen_stl(filename)
        return y0 * scale, y1 * scale

    def gen_scad_header(self, f, y0, y1):
        blade_stl_filename = self.param.name + "_blade.stl"
        f.write(
            textwrap.dedent(
                f"""
            center_hole = {self.param.center_hole * 1000.0};
            hub_diameter = {self.param.hub_radius * 2000};
            hub_height = {self.param.hub_depth * 1000.0};
            n_blades = {self.n_blades};
            y_min = {y0};
            y_max = {y1};
            blade_name = "{blade_stl_filename}";
            """
            )
        )

    def gen_scad(self, filename, y0, y1, ccw=False):
        """Create an OpenSCAD file for the propeller"""
        f = open(filename, "w")
        self.gen_scad_header(f, y0, y1)

        f.write(
            textwrap.dedent(
                """
            module blade() {
                import(blade_name);
            }

            module hub() {
                difference() {
                    cylinder(d=hub_diameter+0.1, h=hub_height, $fn=61);
                    cylinder(d=center_hole, h=55, center=true, $fn=31);
                }
            }
            module prop() {
                union() {
                    for(angle = [0 : (360/n_blades) : 360]) {
                        rotate(angle) blade();
                    }
                    translate([0,0,-hub_height+y_max]) hub();
                }
            }
            """
            )
        )
        if ccw:
            f.write("mirror([1,0,0]) prop();\n")
        else:
            f.write("prop();\n")
        f.close()

    def gen_removable_blade_scad(self, filename, y0, y1, ccw=False):
        """Create an OpenSCAD file for the propeller"""
        f = open(filename, "w")
        self.gen_scad_header(f, y0, y1)

        # template_file = open('blade_template.scad', 'r')
        template_file = pkg_resources.open_text(templates, "blade_template.scad")

        template = template_file.read()
        template_file.close()

        f.write(template)
        f.close()

    def tip_loss(self, r, phi):
        f = (self.n_blades * (self.param.radius - r * 0.96)) / (2.0 * r * np.sin(phi))
        tip_loss = 2.0 * np.arccos(np.exp(-f)) / np.pi

        f = (self.n_blades * (r - self.param.hub_radius * 0.95)) / (
            2.0 * r * np.sin(phi)
        )
        hub_loss = 2.0 * np.arccos(np.exp(-f)) / np.pi
        return tip_loss * hub_loss

    def full_optimize(self, optimum_torque, optimum_rpm, thrust):
        self.blade_elements = []
        u_0 = self.param.forward_airspeed

        dv_goal = optimize.dv_from_thrust(thrust, R=self.param.radius, u_0=u_0)
        radial_points = np.linspace(
            self.param.radius, self.param.hub_radius, self.radial_steps
        )

        total_thrust = 0.0
        total_torque = 0.0
        omega = (optimum_rpm / 60.0) * 2.0 * np.pi
        dr = abs(radial_points[0] - radial_points[1])

        twist_angles = []
        chords = []
        theta = 0.0  # start guess
        dv = dv_goal  # start guess
        a_prime = 0.001  # start guess
        prev_twist = 0.0

        import matplotlib.pyplot as plt

        phi_all = np.arctan((u_0 + dv_goal) / (omega * radial_points))
        plt.plot(radial_points, self.tip_loss(radial_points, 0.1), label="tip loss")
        plt.plot(radial_points, phi_all, label="phi")
        plt.plot(radial_points, np.cos(phi_all), label="cos(phi)")
        plt.legend()
        plt.grid(True)
        plt.xlabel("radius")
        plt.show()
        # return None
        for r in radial_points:
            u = u_0 + dv_goal
            v = omega * r
            phi = np.arctan(u / v)

            dv_modified = dv_goal * self.tip_loss(r, phi)
            be = self.new_foil(r, optimum_rpm, prev_twist)
            # x, fun = optimize.design_for_dv(foil_simulator=be.fs, dv_goal=dv_modified, \
            # rpm = optimum_rpm, B = self.n_blades, r = r, dr=dr, u_0 = u_0)

            y_limit = self.get_max_depth(r)
            x_limit = self.get_max_chord(r, prev_twist)
            maxchord = be.foil.get_max_chord(
                x_limit, y_limit, prev_twist
            )  # Assumes that the foil chord is 1.0

            x, fun = optimize.optimize_all(
                foil_simulator=be.fs,
                dv_goal=dv_modified,
                rpm=optimum_rpm,
                B=self.n_blades,
                r=r,
                dr=dr,
                u_0=u_0,
                maxchord=maxchord,
            )
            theta, dv, a_prime, chord = x
            be.set_chord(chord)
            # if (fun > 0.03):
            # logger.info("Rescan around {}".format(np.degrees(phi)))
            # opt = 9999.9
            # for th_deg in np.arange(np.degrees(phi)-5, np.degrees(phi)+15, 0.5):
            # dv_test, a_prime_test, err = optimize.bem_iterate(foil_simulator=be.fs, \
            # dv_goal=dv_modified, theta = np.radians(th_deg), \
            # rpm = optimum_rpm, B = self.n_blades, r = r, dr=dr, u_0 = u_0)
            # logger.info("err={:5.4f}, th={:4.2f}, dv={:4.2f}, a'={:5.4f}".format(err, th_deg, dv_test, a_prime_test))
            # delta = abs(dv_test - dv_modified)
            # if (err < 0.01) and (delta < opt):
            # opt = delta
            # dv = dv_test
            # a_prime = a_prime_test
            # theta = np.radians(th_deg)
            # logger.info("Best Delta {} at {} deg".format(delta, th_deg))

            ##x, fun = optimize.design_for_dv(foil_simulator=be.fs, \
            ##th_guess=theta, dv_guess=dv, a_prime_guess=a_prime, dv_goal=dv_modified, \
            ##rpm = optimum_rpm, B = self.n_blades, r = r, dr=dr, u_0 = u_0)
            ##if (fun < 0.01):
            ##theta, dv, a_prime = x

            be.set_twist(theta)
            be.set_bem(dv, a_prime)
            twist_angles.append(theta)
            chords.append(be.foil.chord)

            dT = be.dT()
            dM = be.dM()
            total_thrust += dT
            total_torque += dM

            logger.info(
                "r={} theta={}, dv={}, a_prime={}, thrust={}, torque={}, eff={} ".format(
                    r, np.degrees(theta), dv, a_prime, dT, dM, dT / dM
                )
            )
            print(be)

            self.blade_elements.append(be)
            prev_twist = theta

        self.blade_elements.reverse()
        twist_angles.reverse()
        chords.reverse()
        # Now smooth the twist angles
        # Now smooth the optimum angles of attack
        twist_angles = np.array(twist_angles)
        chords = np.array(chords)
        coeff = np.polyfit(radial_points[::-1], twist_angles, 4)
        twist_angle_poly = np.poly1d(coeff)

        # coeff = np.polyfit(radial_points[::-1], chords, 4)
        # chord_poly = np.poly1d(coeff)
        c_points = np.concatenate(
            (
                np.array([0, self.param.hub_radius / 2, 0.9 * self.param.hub_radius]),
                radial_points[::-1],
            )
        )
        extra_chords = np.concatenate(
            (
                0.9
                * np.array(
                    [self.param.hub_depth, self.param.hub_depth, self.param.hub_depth]
                ),
                chords,
            )
        )
        chord_poly = PchipInterpolator(c_points, smooth(extra_chords))

        import matplotlib.pyplot as plt

        plt.plot(radial_points[::-1], np.degrees(twist_angles), label="twist angles")
        plt.plot(
            radial_points[::-1],
            np.degrees(twist_angle_poly(radial_points[::-1])),
            label="Smoothed twist angles",
        )
        plt.plot(c_points, extra_chords * 1000, label="Extra Chords (mm)")
        plt.plot(c_points, chord_poly(c_points) * 1000, label="Chords smoothed")
        plt.legend()
        plt.grid(True)
        plt.xlabel("Radius (m)")
        plt.ylabel("Twist (degrees)")
        # plt.savefig()
        plt.show()
        print("Smoothed Blade Form")

        for be in self.blade_elements:
            be.set_chord(chord_poly(be.r))
            be.set_twist(twist_angle_poly(be.r))
            print(be)

        torque, thrust = self.get_forces(optimum_rpm)
        return torque, thrust


class NACAProp(Prop):
    """Prop that uses NACA Airfoils"""

    def new_foil(self, r, rpm, twist):
        return self.new_blade_element(foil.NACA4, r, rpm, twist)


class ARADProp(Prop):
    """Prop that uses ARAD Airfoils"""

    def new_foil(self, r, rpm, twist):
        from proply import foil_ARA

        return self.new_blade_element(foil_ARA.ARADFoil, r, rpm, twist)
