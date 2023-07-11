import os
import numpy as np
import pandas as pd


# Gets nearest number in array, returns index
# For finding nearest dark-frame exposure
def nearest_number(array, value):
    try:
        array = pd.to_numeric(array)
        value = pd.to_numeric(value)
        array = np.asarray(array)
        index = (np.abs(array - value)).argmin()
        result = array[index]
        return result
    except:
        print('nearest_number error')
        pass

def nearest_date(array, value):
    try:
        return min(array, key=lambda x: abs(x - value))
    except:
        print('nearest_date error')
        pass

def background_remove_list(images_list, dark_frames_list, to_tiff_path):
    try:
        # Get relative path
        images_relative_path = images_list['relative_path']
        n_images = range(0, len(images_list))
        # Add to_tiff_path
        background_remove_image_list = []
        for image_x in n_images:
            image_path = os.path.join(to_tiff_path, images_relative_path[image_x])
            csv_path = os.path.join(image_path, 'channels_metadata.csv')
            if os.path.exists(csv_path):
                channel_metadata = pd.read_csv(csv_path)
                channel_metadata['date'] = images_list['date'][image_x]
                channel_metadata = channel_metadata.assign(
                    image_path=lambda dataframe: dataframe['protein_name'].map(
                        lambda protein_name: os.path.join(image_path, protein_name + '.tif')),
                    exposure=lambda dataframe: dataframe['exposure'].map(
                        lambda exposure: exposure.split(' ')[0] if exposure.split(' ')[1] == 'ms' else exposure)
                )
                background_remove_image_list.append(channel_metadata)
        # Combine rows
        background_remove_image_list = pd.concat(background_remove_image_list)
        # Reset index
        background_remove_image_list = background_remove_image_list.reset_index()
        background_remove_image_list = background_remove_image_list.drop(columns=['index'])
        # Select dark-frame image
        n_protein_images = range(0, len(background_remove_image_list))
        dark_frame_images = []
        for protein_image_x in n_protein_images:
            # Only include images that exist
            if os.path.exists(background_remove_image_list['image_path'][protein_image_x]):
                # Get nearest exposure
                img_df_exposure = nearest_number(dark_frames_list['exposure'],
                                                 background_remove_image_list['exposure'][protein_image_x])
                exposure_dark_frames_list = pd.to_numeric(dark_frames_list['exposure']) == img_df_exposure
                exposure_dark_frames_list = dark_frames_list[exposure_dark_frames_list]
                # Get nearest date
                img_df_date = nearest_date(exposure_dark_frames_list['date'],
                                           background_remove_image_list['date'][protein_image_x])
                date_dark_frames_list = exposure_dark_frames_list['date'] == img_df_date
                date_dark_frames_list = exposure_dark_frames_list[date_dark_frames_list]
                date_dark_frames_list = date_dark_frames_list.reset_index()
                dark_frame_image = date_dark_frames_list['path'][0]
                dark_frame_images.append(dark_frame_image)
        background_remove_image_list['dark_frame'] = dark_frame_images

        return n_protein_images, background_remove_image_list
    except:
        print('background_remove_list error')
        pass

def median_remove_list(background_remove_image_list, exclusion_channels):
    try:
        # Image list
        n_protein_images = range(0, len(background_remove_image_list))
        df_rm_images = []
        for protein_image_x in n_protein_images:
            image_path = background_remove_image_list['image_path'][protein_image_x]
            path = os.path.dirname(image_path)
            file = os.path.basename(image_path)
            if not os.path.splitext(file)[0] in exclusion_channels:
                file = os.path.splitext(file)[0] + '_darkframe_removed.tif'
                df_rm_image = os.path.join(path, file)
                if os.path.exists(df_rm_image):
                    df_rm_images.append(df_rm_image)

        return df_rm_images
    except:
        print('median_remove_list error')
        pass

def combine_images_list(median_remove_image_list, table_ending, merge_ending):
    try:
        data = median_remove_image_list
        to_combine = pd.DataFrame(data, columns=['input'])
        to_combine = to_combine.replace(table_ending, merge_ending, regex=True)

        to_combine = to_combine.assign(
            folder=lambda dataframe: dataframe['input'].map(lambda folder: os.path.dirname(folder)),
            output=lambda dataframe: dataframe['folder'].map(
                lambda folder: os.path.join(folder, 'Combined' + merge_ending))
        )

        return to_combine
    except:
        print('median_remove_list error')
        pass
