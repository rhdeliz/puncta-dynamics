'''
Converts ND2 file to multi-page TIFF
'''

import os, shutil, ray, time, tifffile, csv, warnings, numpy, difflib
from pims_nd2 import ND2_Reader  # Install as pims-nd2
from nd2reader import ND2Reader

def flatten_dictionary(dd, separator='_', prefix=''):
    try:
        return {prefix + separator + k if prefix else k: v
                for kk, vv in dd.items()
                for k, v in flatten_dictionary(vv, separator, kk).items()
                } if isinstance(dd, dict) else {prefix: dd}
    except:
        print('flatten_dictionary error')
        pass

@ray.remote
# images_list
def nd2_to_tiff(image_x_path, image_x_cohort, image_x_channels, to_tiff_path, tiff_compression_level):
    # Get image parameters
    try:
        image_x_name = os.path.basename(image_x_path)
        image_x_name = os.path.splitext(image_x_name)[0]
        save_path = os.path.join(to_tiff_path, image_x_cohort, image_x_name)

        if not os.path.exists(save_path):
            os.makedirs(save_path)
    except:
        print('Image parameters error')
        pass

    try:
        warnings.filterwarnings("ignore")
        img = ND2_Reader(image_x_path)
        warnings.filterwarnings("default")
    except:
        print('Cannot open ND2 image')
        pass

    try:
        # Get channel names
        n_channels = img.metadata['plane_count'] - 1

        channel_names = []
        i = 0
        while i <= n_channels:
            plane_x = 'plane_' + str(i)
            name = img.metadata[plane_x]['name']
            channel_names.append(name)
            i += 1

        channel_numbers = range(0, img.metadata['plane_count'])
        channel_numbers = list(channel_numbers)
        # Get frame count
        n_frames = img.sizes['t'] - 1
    except:
        print('Cannot get channel parameters')
        pass

    try:
        metadata_text = ND2Reader(image_x_path)
        warnings.filterwarnings("ignore")
        timesteps = metadata_text.timesteps
        warnings.filterwarnings("default")
        metadata_text = metadata_text.parser._raw_metadata.image_text_info
        metadata_text = metadata_text[b'SLxImageTextInfo']
        metadata_text = metadata_text[b'TextInfoItem_5']
        metadata_text = metadata_text.decode()

        # Save metadata to text file
        file = os.path.join(save_path, 'metadata.txt')
        file = open(file, 'w')
        file.write(metadata_text)
        file.close()

        # Save actual times
        file = os.path.join(save_path, 'timesteps.csv')
        numpy.savetxt(file, timesteps, delimiter=",")

        # Save metadata to csv
        img_metadata = flatten_dictionary(img.metadata)
        with open(os.path.join(save_path, 'metadata.csv'), 'w') as f:
            f.write('parameter,value,value2,value3\n')
            for key in img_metadata.keys():
                f.write('%s,%s\n' % (key, img_metadata[key]))

        # Get channel emmisions
        emissions = []
        for channel_number in channel_numbers:
            emission = 'plane_' + str(channel_number) + '_emission_nm'
            emission = img_metadata[emission]
            emission = str(int(emission))
            emissions.append(emission)

        # Get select image metadata
        search_terms = ['\r\n  N-STORM Angle: ', '\r\n  N-STORM Direction: ', '\r\n  N-STORM Focus: ', ]
        tirf_parameters = []
        try:
            for term in search_terms:
                tirf_parameter = metadata_text.split(term)[1]
                tirf_parameter = tirf_parameter.split('\r\n')[0]
                tirf_parameter = tirf_parameter.split(' ')[0]
                tirf_parameters.append(tirf_parameter)
        except:
            print('Metadata missing information')
            # Placeholder
            tirf_parameters = 0,0,0
            pass

        # Get illumination data
        exposures = []
        powers = []
        excitations = []
        for channel_name in channel_names:
            try:
                # Get metadata
                separator = '\r\n Name: ' + channel_name
                channel_data = metadata_text.split(separator)[1]
                # Get exposure
                exposure = channel_data.split('\r\n  Exposure: ')[1]
                exposure = exposure.split('\r\n')[0]
                exposures.append(exposure)
                # Get power
                power = channel_data.split('; On; ')[0]
                power = power.split(' Power:')
                power = power[len(power)-1]
                power = power.replace(' ','')
                powers.append(power)
                # Get excitation
                excitation = channel_data.split('; On;')[0]
                excitation = excitation.split('; ExW:')
                excitation = excitation[len(excitation)-1]
                excitation = excitation.replace(' ','')
                excitation = excitation.split(';Power:')[0]
                excitations.append(excitation)
            except:
                print('Channel metadata error')
                pass
        # Get protein names
        protein_names = []
        for channel_name in channel_names:
            try:
                temp_image_x_channels = []
                for c in image_x_channels.index:
                    c = c.upper()
                    temp_image_x_channels.append(c)

                closest_match = difflib.get_close_matches(channel_name.upper(), temp_image_x_channels, cutoff=0.05)[0]
                closest_match = temp_image_x_channels[temp_image_x_channels.index(closest_match)]

                if closest_match.upper() in channel_name.upper():
                    protein_name = image_x_channels[temp_image_x_channels.index(closest_match)]
                else:
                    protein_name = channel_name
            except:
                protein_name = channel_name
                pass
            protein_names.append(protein_name)
        # Save channel metatata
        channel_metadata = os.path.join(save_path, 'channels_metadata.csv')
        with open(channel_metadata, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['protein_name', 'channel', 'power', 'excitation', 'emmision', 'exposure', 'angle', 'direction', 'focus'])
            for c in channel_numbers:
                writer.writerow([protein_names[c], channel_names[c], powers[c], excitations[c], emissions[c], exposures[c], tirf_parameters[0], tirf_parameters[1], tirf_parameters[2]])
    except:
        print('Metadata error')
        pass

    # Split channels
    for channel_number in channel_numbers:
        try:
            frames = range(0, n_frames)
            img_channel = []
            for t in frames:
                orig = img.get_frame_2D(c=channel_number, x=0, y=0, t=t)
                img_channel.append(orig)
            # Save file
            out_name = protein_names[channel_number] + '.tif'
            out_path = os.path.join(save_path, out_name)
            tifffile.imsave(out_path, img_channel, bigtiff=True, compress=tiff_compression_level,
                            dtype=img_channel[0].dtype)
            time.sleep(5)
        except:
            print('Cannot split channels')
            pass

    try:
        # Move ND2 file once finished
        shutil.move(image_x_path, save_path)
        # Close image
        img.close()
        time.sleep(5)
        return image_x_name
    except:
        print('Cannot move processed nd2 file')
        pass
