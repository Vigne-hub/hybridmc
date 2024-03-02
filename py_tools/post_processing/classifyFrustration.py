import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from catboost import CatBoostClassifier
from catboost import Pool
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold

from scipy.stats import randint
from sklearn.metrics import accuracy_score
from sklearn.tree import plot_tree
from numpy import mean
from numpy import std
import matplotlib.pyplot as plt



tree = False
randomForest = False
gradientBoost = False
catboost = True
adaboost = False
# Read the CSV file into a pandas DataFrame
data = pd.read_csv('frustration.csv')

# Split the data into features (X) and target (y)
X = data.iloc[:, :-1]  # Features are all columns except the last one
y = data.iloc[:, -1]   # Target is the last column

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

#print('X is ', X)
# Create a Decision Tree Classifier
if gradientBoost:
    # Define the hyperparameter grid
    param_dist = {
        'n_estimators': randint(100, 200),
        'learning_rate': [0.01, 0.05, 0.1, 0.5, 1],
        'max_depth': [3, 4, 5, 6, 7, 8],
        'min_samples_split': randint(2, 20),
        'min_samples_leaf': randint(1, 30)
    }

    # Create a Gradient Boosting classifier
    gb_classifier = GradientBoostingClassifier(random_state=42)

    # Create RandomizedSearchCV object
    gb_random_search = RandomizedSearchCV(gb_classifier, param_distributions=param_dist, n_iter=100, cv=3,
                                          random_state=42, n_jobs=-1)

    # Fit RandomizedSearchCV to training data
    gb_random_search.fit(X_train, y_train)

    # Get the best parameters
    best_params = gb_random_search.best_params_
    print("Best Boost Parameters:", best_params)

    # Predict on the test set
    y_pred = gb_random_search.predict(X_test)
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Boost Accuracy:", accuracy)

    # Get feature importances from the trained model
    feature_importances = gb_random_search.best_estimator_.feature_importances_

    # Sort feature importances in descending order
    sorted_indices = np.argsort(feature_importances)[::-1]

    # Extract feature names
    feature_names = [f'Feature {i + 1}' for i in range(X.shape[1])]

    # Plot feature importances
    plt.figure(figsize=(10, 6))
    plt.bar(range(X.shape[1]), feature_importances[sorted_indices], align='center')
    plt.xticks(range(X.shape[1]), np.array(feature_names)[sorted_indices], rotation=90)
    plt.xlabel('Features')
    plt.ylabel('Feature Importance')
    plt.title('Feature Importances')
    plt.savefig('boost.png')


if randomForest:
    # Define the hyperparameter grid
    param_dist = {
        'n_estimators': randint(50, 120),
        'max_depth': [None] + list(randint(3, 20).rvs(10)),
        #'max_features': ['auto', 'sqrt', None] + list(randint(1, 20).rvs(10)),
        'min_samples_split': randint(2, 20),
        'min_samples_leaf': randint(1, 20),
        'bootstrap': [True, False]
    }

    # Create a Random Forest classifier
    rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

    # Create RandomizedSearchCV object
    rf_random_search = RandomizedSearchCV(rf_classifier, param_distributions=param_dist, n_iter=100, cv=3,
                                          random_state=42, n_jobs=-1)

    # Fit RandomizedSearchCV to training data
    rf_random_search.fit(X_train, y_train)

    # Get the best parameters
    best_params = rf_random_search.best_params_
    print("Best Parameters:", best_params)

    # Predict on the test set
    y_pred = rf_random_search.predict(X_test)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy:", accuracy)


    # Train the classifier
    rf_classifier.fit(X_train, y_train)

    # Predict on the test set
    y_pred = rf_classifier.predict(X_test)

    print("Random forest Accuracy:", accuracy)

    # Get feature importances from the trained model
    feature_importances = rf_random_search.best_estimator_.feature_importances_

    # Sort feature importances in descending order
    sorted_indices = np.argsort(feature_importances)[::-1]

    # Extract feature names
    feature_names = [f'Feature {i + 1}' for i in range(X.shape[1])]

    # Plot feature importances
    plt.figure(figsize=(10, 6))
    plt.bar(range(X.shape[1]), feature_importances[sorted_indices], align='center')
    plt.xticks(range(X.shape[1]), np.array(feature_names)[sorted_indices], rotation=90)
    plt.xlabel('Features')
    plt.ylabel('Feature Importance')
    plt.title('Feature Importances')
    plt.savefig('forest.png')

if tree:
    clf = DecisionTreeClassifier(max_leaf_nodes=20, random_state=0)
    clf.fit(X_train, y_train)

    # Predict the target values for the test set
    y_pred = clf.predict(X_test)
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Decision tree Accuracy:", accuracy)

if adaboost:
    # Create a weak learner (Decision Tree Classifier)
    weak_learner = DecisionTreeClassifier(max_depth=1)  # Stump
    # Create AdaBoostClassifier with the weak learner
    clf = AdaBoostClassifier(base_estimator=weak_learner, n_estimators=50, random_state=42)
    # Train the classifier on the training data
    clf.fit(X_train, y_train)

    # Predict the target values for the test set
    y_pred = clf.predict(X_test)

    # Calculate the accuracy of the classifier
    accuracy = accuracy_score(y_test, y_pred)
    print("Ada Accuracy:", accuracy)

if tree:
    n_nodes = clf.tree_.node_count
    children_left = clf.tree_.children_left
    children_right = clf.tree_.children_right
    feature = clf.tree_.feature
    threshold = clf.tree_.threshold
    values = clf.tree_.value

    node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
    is_leaves = np.zeros(shape=n_nodes, dtype=bool)
    stack = [(0, 0)]  # start with the root node id (0) and its depth (0)
    while len(stack) > 0:
        # `pop` ensures each node is only visited once
        node_id, depth = stack.pop()
        node_depth[node_id] = depth

        # If the left and right child of a node is not the same we have a split
        # node
        is_split_node = children_left[node_id] != children_right[node_id]
        # If a split node, append left and right children and depth to `stack`
        # so we can loop through them
        if is_split_node:
            stack.append((children_left[node_id], depth + 1))
            stack.append((children_right[node_id], depth + 1))
        else:
            is_leaves[node_id] = True

    print(
        "The binary tree structure has {n} nodes and has "
        "the following tree structure:\n".format(n=n_nodes)
    )
    for i in range(n_nodes):
        if is_leaves[i]:
            print(
                "{space}node={node} is a leaf node with value={value}.".format(
                    space=node_depth[i] * "\t", node=i, value=values[i]
                )
            )
        else:
            print(
                "{space}node={node} is a split node with value={value}: "
                "go to node {left} if X[:, {feature}] <= {threshold} "
                "else to node {right}.".format(
                    space=node_depth[i] * "\t",
                    node=i,
                    left=children_left[i],
                    feature=feature[i],
                    threshold=threshold[i],
                    right=children_right[i],
                    value=values[i],
                )
            )
    plot_tree(clf)
    plt.savefig('decision_tree.png')
    #plt.show()


if catboost:
    # Create a CatBoost Pool object
    train_pool = Pool(X_train, label=y_train)

    catboost_classifier = CatBoostClassifier(iterations=1000,  # Number of trees (boosting iterations)
                                             learning_rate=0.1,  # Learning rate
                                             depth=8,  # Depth of the trees
                                             loss_function='Logloss',  # Loss function for binary classification
                                             random_seed=42,  # Random seed for reproducibility
                                             early_stopping_rounds=20,
                                             verbose=10)  # Print progress every 100 iterations

    # Train the classifier
    #catboost_classifier.fit(X_train, y_train, eval_set=(X_test, y_test), plot=True)
    catboost_classifier.fit(train_pool)

    # Predict on the test set
    y_pred = catboost_classifier.predict(X_test)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Catboot Accuracy:", accuracy)

    # Get feature importance scores
    feature_importance = catboost_classifier.get_feature_importance(train_pool)
    sorted_indices = np.argsort(feature_importance)[::-1]

    # Print sorted feature importance scores
    print("Sorted Feature Importance Scores:")
    for i, feature_index in enumerate(sorted_indices):
        print(f"Feature {feature_index + 1}: {feature_importance[feature_index]}")



