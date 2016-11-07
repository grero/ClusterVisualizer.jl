using ClusterVisualizer

"""
Create `nclasses` clusters that evolve in time over `nbins` time bins according to a random walk process.

	function prepare_data(nbins, nn,nclasses=3)
"""
function prepare_data(nbins, nn,nclasses=3)
  X = zeros(Float64, 3,nn, nbins)
	μ = zeros(3,nbins,nclasses)
  μ[:,1,:] = 0.0
	for i in 1:nclasses
		for j in 2:nbins
			μ[:,j,i] = μ[:,j-1,i] + 0.1*randn(3)
		end
	end
	labels = zeros(Int64,nn)
	for i in 1:nn	
		l = rand(1:nclasses)
		labels[i] = l
		for j in 1:nbins
			X[:,i,j] = μ[:,j,l] + 0.5randn(3)
		end
	end
  X,labels
end

X,labels = prepare_data(500,100,3)
ClusterVisualizer.animate_clusters(X,labels,10;show_paths=true)


