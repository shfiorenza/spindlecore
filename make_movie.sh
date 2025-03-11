#usage: ./make_movie.sh <FPS> <sim frames to use per mov frame> <dir/name of sim>
#exmpl: ./make_movie.sh 40 40 frames/test
# SF note: I don't fully understand the difference between inputs 1 and 2
# SF note: however, it seems to work best when both are equal
#some ffmpeg notes:
#"ffmpeg "
        # frame rate (Hz)
        #   "-r 20 "
        # frame size (width x height)
        #    "-s 1080x720 "
        # input files
        #    "-i "
        # image_directory
        #    f"/fig_%04d.{image_type} "
        # video codec
        #    "-vcodec libx264 "
        # video quality, lower means better
        #    "-crf 25 "
        # pixel format
        #    "-pix_fmt yuv420p "
        # gaussian_output file
        #   movie_directory
ffmpeg -y -f image2 -framerate $1 -i $3_%*.bmp -vcodec libx264 -profile baseline -pix_fmt yuv420p -r $2 -q:v 0.8 $3.mov