class UnetParam(object):

    def __init__(self, cnf, unet_img_size=(512, 512)):
        self.from_top_nucleus_unet_model = cnf.from_top_nucleus_unet_model
        self.nucleus_unet_model = cnf.nucleus_unet_model
        self.unet_model_scale = cnf.unet_model_scale
        self.unet_model_thrh = cnf.unet_model_thrh
        self.unet_img_size = unet_img_size

class ImgResolution(object):

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class CellsImg(object):

    def __init__(self, img_name, resolution, cells):
        self.name = img_name
        self.cells = cells
        self.resolution = resolution

class TestStructure(object):

    def __init__ (self, fibers, nodes, pairs, resolution):
        self.fibers = fibers
        self.nodes = nodes
        self.pairs = pairs
        self.resolution = resolution