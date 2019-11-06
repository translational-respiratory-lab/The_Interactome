import pandas
import numpy
def _silhouette_samples(arr, labels):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.utils import check_X_y

    def check_number_of_labels(n_labels, n_samples):
        if not 1 < n_labels < n_samples:
            raise ValueError("Number of labels is %d. Valid values are 2 "
                             "to n_samples - 1 (inclusive)" % n_labels)

    arr, labels = check_X_y(arr, labels, accept_sparse=['csc', 'csr'])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    check_number_of_labels(len(le.classes_), arr.shape[0])

    unique_labels = le.classes_
    n_samples_per_label = numpy.bincount(labels, minlength=len(unique_labels))

    # For sample i, store the mean distance of the cluster to which
    # it belongs in intra_clust_dists[i]
    intra_clust_aff = numpy.zeros(arr.shape[0], dtype=arr.dtype)

    # For sample i, store the mean distance of the second closest
    # cluster in inter_clust_dists[i]
    inter_clust_aff = intra_clust_aff.copy()
    
    
    for curr_label in range(len(unique_labels)):

        # Find inter_clust_dist for all samples belonging to the same
        # label.
        mask = labels == curr_label
        current_distances = arr[mask]
               
        # Leave out current sample.
        n_samples_curr_lab = n_samples_per_label[curr_label] - 1
        
        if n_samples_curr_lab != 0:
            intra_clust_aff[mask] = numpy.sum(
                current_distances[:, mask], axis=1) / n_samples_curr_lab

        # Now iterate over all other labels, finding the mean
        # cluster distance that is closest to every sample.
        for other_label in range(len(unique_labels)):
            if other_label != curr_label:
                other_mask = labels == other_label
                other_distances = numpy.mean(
                    current_distances[:, other_mask], axis=1)
                inter_clust_aff[mask] = numpy.maximum(
                    inter_clust_aff[mask], other_distances)

    sil_samples = intra_clust_aff - inter_clust_aff
    sil_samples /= numpy.maximum(intra_clust_aff, inter_clust_aff)

    # score 0 for clusters of size 1, according to the paper
    sil_samples[n_samples_per_label.take(labels) == 1] = 0

    return sil_samples


def silhouette_score(arr, labels):
	z=_silhouette_samples(arr, labels)
	df=pandas.DataFrame([z,labels])
	df=df.transpose()
	x=numpy.mean(z)
	s_c=[x]
	for i in set(labels):
		c=df[df[1]==i][0]
		s_c.append(numpy.mean(c))
	return s_c
