# Time to draw
from PIL import Image,ImageDraw,ImageFont
import PIL
image = Image.new( 'RGBA', [200,200], ( 0, 0, 0, 0 ) ) # Create a blank picture
drawer = ImageDraw.Draw( image ) # Create a drawer to the picture
drawer.rectangle( [200,10,0,20] ,fill=(100,149,237) )
drawer.rectangle( [20,   #width
                   10,
                   30,    #position horizontal
                   20] ,fill=(167,191,66) )
#font = ImageFont.truetype('/arial.ttf',12)
#drawer.text([200,10,0] ,'GENE', fill=(0,0,0))
drawer.text((55, 0),'sequencia teste',fill=(0,0,0))

image.save( 'teste.png', 'PNG' )


it = box_step # * 7
for p in data:

    drawer.text( ( p[ 'start' ], it ), p[ 'name' ], fill=outline_color ) # Draw the name of the feature

    box = ( p[ 'start' ], it + font_size[ 1 ], # Feature's box
    p[ 'end' ], it + box_heigh + font_size[ 1 ] )

    drawer.rectangle( box, fill=outline_color ) # Out box
    drawer.rectangle( ( box[ 0 ] + border_width, box[ 1 ] + border_width, box[ 2 ] - border_width, box[ 3 ] - border_width ),fill=fill_colors[ p[ 'type' ] ] ) # In box

    it += box_step + box_heigh + font_size[ 1 ]


# //