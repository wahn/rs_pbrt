import math
# Blender
import bpy
import bpy_extras
import mathutils

# define function
def pbrt_sphere(radius = 1.0, zmin = -1.0, zmax = 1.0, phimax = 360.0, steps = 16):
    # default values?
    if radius == 1.0 and zmin == -1.0 and zmax == 1.0 and phimax == 360.0:
        print("pbrt_sphere(radius = %s, zmin = %s, zmax = %s, phimax = %s)" %
              (radius, zmin, zmax, phimax))
        # create ico sphere
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=4)
        # shading: smooth
        bpy.ops.object.shade_smooth()
        # change name(s)
        obj = bpy.context.selected_objects[0]
        obj.name = "PbrtSphere"
        obj.data.name = "PbrtSphere"
        # custom properties
        obj = bpy.data.objects["PbrtSphere"]
        obj['_RNA_UI'] = {}
        obj["radius"] = radius
        obj["zmin"] = zmin
        obj["zmax"] = zmax
        obj["phimax"] = phimax
        return
    elif phimax == 360.0:
        # proper sphere?
        if zmin == -radius and zmax == radius:
            print("pbrt_sphere(radius = %s, zmin = %s, zmax = %s, phimax = %s)" %
                  (radius, zmin, zmax, phimax))
            # create ico sphere (with radius)
            bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=4, size=radius)
            # shading: smooth
            bpy.ops.object.shade_smooth()
            # change name(s)
            obj = bpy.context.selected_objects[0]
            obj.name = "PbrtSphere"
            obj.data.name = "PbrtSphere"
            # custom properties
            obj = bpy.data.objects["PbrtSphere"]
            obj['_RNA_UI'] = {}
            obj["radius"] = radius
            obj["zmin"] = zmin
            obj["zmax"] = zmax
            obj["phimax"] = phimax
            return
    # check arguments (and correct them)
    if radius <= 0.0:
        print("WARNING: radius <= 0.0 (%s <= 0.0)" % radius)
        new_radius = 1.0
    else:
        new_radius = radius
    if zmin < -new_radius:
        print("WARNING: zmin < -radius (%s < %s)" % (zmin, -new_radius))
        new_zmin = -new_radius
    else:
        new_zmin = zmin
    if zmax > new_radius:
        print("WARNING: zmax > radius (%s > %s)" % (zmax, new_radius))
        new_zmax = new_radius
    else:
        new_zmax = zmax
    if phimax <= 0.0:
        print("WARNING: phimax <= 0.0 (%s <= 0.0)" % phimax)
        new_phimax = 360.0
    elif phimax > 360.0:
        print("WARNING: phimax > 360.0 (%s > 360.0)" % phimax)
        new_phimax = 360.0
    else:
        new_phimax = phimax
    print("pbrt_sphere(radius = %s, zmin = %s, zmax = %s, phimax = %s)" %
          (new_radius, new_zmin, new_zmax, new_phimax))
    # some edges to spin
    verts = []
    edges = []
    start = math.degrees(math.acos(-new_zmin/new_radius))
    stop = math.degrees(math.acos(-new_zmax/new_radius))
    step = (stop - start) / float(steps)
    count = 0
    angle = start
    while angle <= stop:
        x = math.sin(math.radians(angle)) * new_radius
        y = 0
        z = -math.cos(math.radians(angle)) * new_radius
        verts.append(mathutils.Vector((x, y, z)))
        if count > 0:
            edges.append([count-1, count])
        count += 1
        angle += step
    faces = []
    mesh = bpy.data.meshes.new(name="PbrtSphere")
    mesh.from_pydata(verts, edges, faces)
    # mesh.validate(verbose=True)
    bpy_extras.object_utils.object_data_add(bpy.context, mesh)
    # spin
    bpy.ops.object.editmode_toggle()
    bpy.ops.mesh.spin(steps=2*steps, angle=-math.radians(new_phimax), axis=(0,0,1))
    bpy.ops.object.editmode_toggle()
    # shading: smooth
    bpy.ops.object.shade_smooth()
    # custom properties
    obj = bpy.data.objects["PbrtSphere"]
    obj['_RNA_UI'] = {}
    obj["radius"] = new_radius
    obj["zmin"] = new_zmin
    obj["zmax"] = new_zmax
    obj["phimax"] = new_phimax

# call function
radius = 1.0
zmin = -0.9
zmax = 0.8
phimax = 270.0
steps = 16
pbrt_sphere(radius = radius, zmin = zmin, zmax = zmax, phimax = phimax, steps = steps)
