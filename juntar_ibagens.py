import sys
from PIL import Image

path = "toy/img/"

ibages_juntadas = 0


def branco():
    new_im = Image.new('RGB', (640, 480), (255, 255, 255))
    new_im.save('toy/img/master/branco.jpg')
    new_im.save('toy/img/dual/branco.jpg')


def juntar_ibages(codigo_master, codigo_dual):
    if type(codigo_master) == int:
        codigo_master = str(codigo_master).zfill(4) + '.png'
    else:
        codigo_dual += '.jpg'
    if type(codigo_dual) == int:
        codigo_dual = str(codigo_dual).zfill(4) + '.png'
    else:
        codigo_dual += '.jpg'
    global ibages_juntadas
    images = [Image.open(path + x) for x in ['master/' + codigo_master,
                                             'dual/' + codigo_dual]]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    ibages_juntadas += 1
    new_im.save('toy/img/ibages/ibagem_' + str(ibages_juntadas).zfill(4) + '.jpg')


if __name__ == '__main__':
    codigos = [(i, 'branco') for i in range(44)]
    codigos += [(43, j) for j in range(44, 58)]
    codigos += [(1058, 1059),
                (1060, 1061),
                (2062, 2063),
                (2064, 2065),
                (2066, 2067),
                (2068, 2069),
                (2070, 2071),
                (2072, 2073),
                (3074, 3075),
                (3076, 3078),
                (3077, 3079),
                (3080, 3081),
                (3082, 3083),
                (3084, 3086),
                (3085, 3087),
                (3088, 3089),
                (3090, 3092),
                (3091, 3093),
                (3094, 3095),
                (3096, 3098),
                (3097, 3099),
                (3100, 3101),
                (3102, 3103),
                (3104, 3105),
                (3106, 3107),
                (3108, 3111),
                (3109, 3112),
                (3110, 3113),
                (4114, 4118),
                (4115, 4119),
                (4116, 4120),
                (4117, 4121),
                (4117, 5122),
                (6123, 6124),
                (6125, 6126),
                ]
    codigos += [(x, x+1) for x in range(6127, 6150, 2)]

    for c in codigos:
        juntar_ibages(*c)
