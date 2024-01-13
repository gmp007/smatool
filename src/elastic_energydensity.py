"""
  SMATool -- Automated toolkit for computing zero and finite-temperature strength of materials

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  Email: cekuma1@gmail.com

""" 
import os
import pandas as pd
import numpy as np
import joblib
from ase import Atoms
import re
import periodictable
from sklearn.model_selection import train_test_split, ShuffleSplit, cross_val_score
#from sklearn.ensemble import RandomForestRegressor, DecisionTreeRegressor, ExtraTreesRegressor, GradientBoostingRegressor, AdaBoostRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
#from xgboost import XGBRFRegressor
from catboost import CatBoostRegressor
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from pymatgen.core.composition import Composition
from matminer.featurizers.composition import ElectronAffinity, ElementFraction
from read_write import simplify_formula
from numpy import cross, linalg



def predict_thickness_2D(atoms, dir_modeldsave):
    # Load the existing data
    existing_data = load_thickness()

    # Get the chemical formula of the new data point
    chem_formula = atoms.get_chemical_formula(mode='hill', empirical=False)
    chem_formula = simplify_formula(chem_formula)

    # Check if the material already exists in the existing data
    if chem_formula in existing_data['MaterialName'].values:
        existing_thickness = existing_data[existing_data['MaterialName'] == chem_formula]['Thickness_Ang'].iloc[0]
        print(f"Thickness for {chem_formula} already exists. Using existing thickness: {existing_thickness} Å.")
        return existing_thickness

    
    # Process and scale the existing data for training
    processed_existing_data = process_dataframe(existing_data)
    X_scaled_existing, y_existing = scale_dataframe(processed_existing_data)
    train_features = X_scaled_existing.columns.tolist()

    # Process and scale the single data point for prediction
    processed_single_data = process_dataframe(pd.DataFrame([chem_formula], columns=['MaterialName']))
    X_scaled_single, _ = scale_dataframe(processed_single_data)

    for feature in train_features:
        if feature not in X_scaled_single.columns:
            X_scaled_single[feature] = 0  # Add missing feature with a value of 0
    X_scaled_single = X_scaled_single[train_features]

    # Load saved model if available, else train a new one
    saved_model_path = os.path.join(dir_modeldsave, 'best_thickness_model.pkl')
    if os.path.exists(saved_model_path):
        model = joblib.load(saved_model_path)
        print("Loaded saved model for predictions.")
    else:
        # Train using the existing bulk data
        model_metrics, model = train_and_save_best_model(X_scaled_existing, y_existing, dir_modeldsave)
        print("Trained a new model for predictions.\n ", model_metrics)

    # Make predictions on the new single data point
    predictions = model.predict(X_scaled_single)

    # Return the predictions
    return predictions




def calculate_energy_densities(max_strain, max_stress, atoms, dim, h_Angstrom):
    # Constants
    eV_to_Joules = 1.602E-19
    m3_to_liters = 1000
    J_to_MJ = 1E-6
    amu_to_kg = 1.66053906660e-27
    Na = 6.02E+23

    # Calculate volume
    if dim == '2D':
        cell = atoms.get_cell()
        area =linalg.norm(cross(cell[0], cell[1]))
        volume = area * h_Angstrom
        strain_energy_density = 0.5*max_strain*max_stress/h_Angstrom*1E10 #(Note stress in 2D is N/m)
    elif dim in ['3D', '1D']:
        volume = atoms.get_volume()
        strain_energy_density = 0.5*max_strain*max_stress*1E9
    volume_m3 = volume * 1e-30  # Convert to m³
    E_eV = strain_energy_density/eV_to_Joules*volume_m3
    density_kg_m3 = sum(atoms.get_masses())* amu_to_kg/volume_m3 
    # Energy calculations
    number_atoms = atoms.get_global_number_of_atoms()
    E_eV_per_atom = E_eV / number_atoms

    E_Joules = E_eV_per_atom * eV_to_Joules

    # Calculate energy density in Joules per cubic meter
    energy_density_J_m3 = E_Joules / volume_m3
    # Convert energy density to MJ/L
    energy_density_MJ_L = energy_density_J_m3 * J_to_MJ / m3_to_liters

    # Convert J/m³ to Wh/m³
    energy_density_Wh_m3 = energy_density_J_m3 / 3600
    
    # Convert Wh/m³ to Wh/kg
    energy_density_Wh_kg = energy_density_Wh_m3 / density_kg_m3

    return energy_density_MJ_L, energy_density_Wh_kg
    


def process_dataframe(data):
    data = pd.DataFrame(data)
    # Create Composition objects
    Comp = [Composition(value) for value in data["MaterialName"]]
    data.loc[:, 'Composition'] = Comp

    # Featurize the dataframe using ElementFraction
    ef = ElementFraction()
    data = ef.featurize_dataframe(data, 'Composition')
    data = data.loc[:, (data != 0).any(axis=0)]

    # Define the molecular_weight function
    def molecular_weight(formula):
        try:
            elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
            weight = 0
            for element, count in elements:
                count = int(count) if count else 1
                weight += getattr(periodictable, element).mass * count
            return weight
        except AttributeError:
            print(f"Warning: Element not found in formula {formula}")
            return None

    # Calculate molecular weights and add them as a new column
    data['MolWeight'] = data["MaterialName"].apply(molecular_weight)

    return data
    


def scale_dataframe(df):
    scaler = StandardScaler()

    # Check if 'Thickness_Ang' is in the DataFrame
    if 'Thickness_Ang' in df.columns:
        # Separate features and target when 'Thickness_Ang' is present
        X = df.drop(columns=["Thickness_Ang", "MaterialName", "Composition"], axis=1)
        y = df["Thickness_Ang"]
    else:
        # If 'Thickness_Ang' is not in the DataFrame, consider all columns as features
        X = df.drop(columns=["MaterialName", "Composition"], axis=1)
        y = None  # y is not available

    # Scale the features
    X_scaled = scaler.fit_transform(X)
    X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

    return X_scaled_df, y



def train_and_save_best_model(X, Y, directory, test_size=0.25, rndseem=101):

    X_train_n, X_test, y_train_n, y_test = train_test_split(X, Y, test_size=test_size, random_state=rndseem)

    augmenter = DataAugmenter(mean=0, std=0.1, num_augmented_samples=10)
    X_train, y_train = augmenter.augment_and_shuffle(X_train_n, y_train_n)


    # Define models and their hyperparameters
    hyper_params_xgb = {
        'objective': 'reg:squarederror',
        'n_estimators': 1000,
        'eval_metric': mean_absolute_error
    }

    MLA = [
      #  RandomForestRegressor(),
      #  DecisionTreeRegressor(),
      #  ExtraTreesRegressor(),
      #  XGBRFRegressor(**hyper_params_xgb),
      #  AdaBoostRegressor(),
      #  GradientBoostingRegressor(),
      #  ExtraTreesRegressor(),
        CatBoostRegressor(loss_function='RMSE', silent=True)
    ]

    # Initialize variables to store results
    algorithms = pd.DataFrame()
    best_model = None
    best_score = -float('inf')

    cv = ShuffleSplit(n_splits=10, test_size=test_size + 0.1, random_state=rndseem)

    for model in MLA:
        try:
            Alg = model.__class__.__name__
            if Alg == 'XGBRFRegressor':
                model.fit(X_train, y_train, eval_set=[(X_test, y_test)])
            else:
                model.fit(X_train, y_train)

            cross_validation = cross_val_score(model, X_train, y_train, cv=cv, scoring='r2')
            mean_cv_score = cross_validation.mean()
            pred = model.predict(X_test)
            adj_R2 = 1 - (1 - model.score(X_train, y_train)) * (len(y_train) - 1) / (len(y_train) - X_train.shape[1] - 1)
            Train_Score = model.score(X_train, y_train)
            Test_Score = model.score(X_test, y_test)
            MSE = metrics.mean_squared_error(y_test, pred)
            MAE = metrics.mean_absolute_error(y_test, pred)
            STD = cross_validation.std()

            if mean_cv_score > best_score:
                best_score = mean_cv_score
                best_model = model

            algorithms = algorithms.append({
                'Algorithm': Alg,
                'Model-Sc': round(Train_Score * 100, 2),
                'Adj-Sc': round(adj_R2 * 100, 2),
                'Test-Sc': round(Test_Score * 100, 2),
                'CV-Sc': round(mean_cv_score * 100, 2),
                'MSE': round(MSE, 2),
                'MAE': round(MAE, 2),
                'STD': round(STD, 2)
            }, ignore_index=True)

        except Exception as e:
            print(f"Exception occurred in {Alg}: {e}")
            pass

    # Save the best model
    if best_model is not None:
        joblib.dump(best_model, f'{directory}/best_thickness_model.pkl')

    return algorithms , best_model

def load_thickness():

    data = {
        "MaterialName": ["C", "MoS2", "WS2", "MoSe2", "WSe2", "BN", "MoTe2", "WTe2", "NbSe2", "TaS2", 
                          "TaSe2", "Bi2Se3", "Bi2Te3", "SnS2", "SnSe2", "HfS2", "HfSe2", "ZrS2", "ZrSe2", "TiS2", 
                          "TiSe2", "ReS2", "ReSe2", "CuI", "P", "Bi2Te3", "NbSe2", "SnS2", "MoTe2", "WTe2", 
                          "NbSe2", "TaS2", "TaSe2", "FeSe", "NiTe2", "PtS2", "PtSe2", "PtTe2", "PdTe2", "VS2", 
                          "VSe2", "MnSe2", "CrS2", "CrSe2", "FePS3", "WS2", "InSe", "SiC", "Si", "Ge",
                          "ZnO", "ZnS", "ZnSe", "CdO", "CdS", "CdSe", "MgO",
                          "PbTe", "PbSe", "SnSe", "SnTe", "GeTe", "SiTe"],
        "Thickness_Ang": [3.43, 6.5, 6.2, 6.7, 6.4, 3.4, 10.8, 9.9, 6.7, 6.4, 
                          6.6, 15, 16, 6.1, 6.5, 6.5, 6.8, 6.4, 6.7, 6, 
                          6.3, 7.2, 7, 2.8, 6.5, 12, 8, 5.8, 10.8, 9.9, 
                          6.7, 6.4, 6.6, 6.2, 7.27, 4.9, 5.9, 6.8, 6.8, 5.9, 
                          6.32, 6.8, 6.16, 6.6, 7.9, 6.2, 6.1, 3.1, 0.67, 0.67,
                          3.05, 3.45, 3.75, 3.4, 3.75, 3.95, 2.3,
                          6.35,  # Average of 6.2 - 6.5
                          5.95,  # Average of 5.8 - 6.1
                          5.55,  # Average of 5.4 - 5.7
                          5.85,  # Average of 5.7 - 6.0
                          5.15,  # Average of 5.0 - 5.3
                          4.85   # Average of 4.7 - 5.0                          
                          ]
        }



    data = pd.DataFrame(data)
    df_avg = data.groupby('MaterialName', as_index=False)['Thickness_Ang'].mean()
    data = df_avg.drop_duplicates(subset=['MaterialName'])

    #data = data.set_index('MaterialName')

    return data


# Augment data to sparse dataset to improve performance
class DataAugmenter:
    def __init__(self, mean=0, std=0.10, num_augmented_samples=2):
        self.mean = mean
        self.std = std
        self.num_augmented_samples = num_augmented_samples

    def add_gaussian_noise(self, data):
        """Add Gaussian noise to data."""
        noise = np.random.normal(self.mean, self.std, data.shape)
        return data + noise

    def augment_data_continuous(self, X, y):
        """Generate noisy versions of each sample for continuous data."""
        augmented_X = []
        augmented_y = []

        for idx, (sample, label) in enumerate(zip(X.values, y)):
            for _ in range(self.num_augmented_samples):
                noisy_sample = self.add_gaussian_noise(sample)
                augmented_X.append(noisy_sample)
                augmented_y.append(label)

        # Create DataFrame from augmented data
        augmented_X_df = pd.DataFrame(augmented_X, columns=X.columns)
        if isinstance(y, pd.Series):
            augmented_y_df = pd.Series(augmented_y, name=y.name)
        else:
            augmented_y_df = pd.DataFrame(augmented_y, columns=y.columns)

        return augmented_X_df, augmented_y_df

    def shuffle_dataset(self, X, y):
        """Shuffle data and labels."""
        combined = pd.concat([X, y], axis=1)
        shuffled = combined.sample(frac=1).reset_index(drop=True)
        return shuffled[X.columns], shuffled[y.name] if isinstance(y, pd.Series) else shuffled[y.columns]

    def augment_and_shuffle(self, X, y):
        augmented_X, augmented_y = self.augment_data_continuous(X, y)
        return self.shuffle_dataset(augmented_X, augmented_y)
