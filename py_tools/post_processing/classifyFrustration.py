import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import plot_tree
import matplotlib.pyplot as plt


tree = True
#tree = False
# Read the CSV file into a pandas DataFrame
data = pd.read_csv('frustration.csv')

# Split the data into features (X) and target (y)
X = data.iloc[:, :-1]  # Features are all columns except the last one
y = data.iloc[:, -1]   # Target is the last column

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)

#print('X is ', X)
# Create a Decision Tree Classifier

if tree:
    clf = DecisionTreeClassifier(max_leaf_nodes=9, random_state=0)
else:
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
print("Accuracy:", accuracy)

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
    plt.show()
else:
    # Get the first weak learner from the AdaBoost ensemble
    weak_learner = clf.estimators_[0]  # Assuming the first estimator is a decision stump

    plt.figure(figsize=(10, 6))
    plot_tree(weak_learner, filled=True, feature_names=X.columns, class_names=True)
    plt.savefig('decision_tree_ada.png')
    plt.show()
