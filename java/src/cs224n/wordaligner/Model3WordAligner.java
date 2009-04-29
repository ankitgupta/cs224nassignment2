/**
 * 
 */
package cs224n.wordaligner;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;


import cs224n.util.Alignment;
import cs224n.util.Counter;
import cs224n.util.CounterMap;
import cs224n.util.Pair;
import cs224n.util.PriorityQueue;
import cs224n.util.SentencePair;

/**
 * @author ankit
 *
 */
public class Model3WordAligner extends WordAligner {

	CounterMap<String, String> tCountMap; //target, source
	CounterMap<String, Integer> nCountMap; //target, phi
	CounterMap<String, String> tProbabilityMap; //source, target
	CounterMap<Integer, String> nProbabilityMap; //phi, target
	private int nbuckets = 6;
	private double[] d, dupdate;
	HashMap<Pair<List<String>, List<String>>, Counter<Alignment>> possibleAlignments;
	CounterMap<List<String>, List<String>> alignmentProbabilities; //source, target
	double p;
	double dmax;
	int phimax;

	HashMap<Integer, List<List<Integer>>> Partitions;
	HashMap<Integer, Integer> Factorial;
	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#alignSentencePair(cs224n.util.SentencePair)
	 */
	@Override
	public Alignment alignSentencePair(SentencePair sentencePair) {
		double probability, max=0;
		Alignment best = null;
		List<String> targetSentence = sentencePair.getEnglishWords();
		List<String> sourceSentence = sentencePair.getFrenchWords();
		//		System.out.println("French "+sourceSentence+" English "+targetSentence);
		for(Alignment alignment : possibleAlignments.get(new Pair<List<String>, List<String>>(sourceSentence, targetSentence)).keySet()) {
			probability = getAlignmentProb(targetSentence, sourceSentence, alignment);
			//			System.out.println(alignment+" "+probability);
			if(max < probability || best==null) {
				probability = max;
				best = alignment;
			}
		}
		//		System.out.println(best);
		return best;
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#getAlignmentProb(java.util.List, java.util.List, cs224n.util.Alignment)
	 */
	@Override
	public double getAlignmentProb(List<String> targetSentence,
			List<String> sourceSentence, Alignment alignment) {
		//		System.out.println(alignmentProbabilities.getCount(sourceSentence, targetSentence));
		return getAlignmentProbHelper(targetSentence, sourceSentence, alignment)/alignmentProbabilities.getCount(sourceSentence, targetSentence);
		//		System.out.println("*******WTF");
		//		return 0;
	}

	private double getAlignmentProbHelper(List<String> targetSentence,
			List<String> sourceSentence, Alignment alignment) {
		int I = targetSentence.size();
		int J = sourceSentence.size();
		double probability = 0.0;
		int phi0 = 0;
		for(int j=0;j<J;j++) {
			int i = alignment.getAlignedTarget(j);
			if(i>0) {
				probability += Math.log(tProbabilityMap.getCount(sourceSentence.get(j), targetSentence.get(i))) + Math.log(d[bucket(i-j*I/J)]);
				// System.out.println(sourceSentence.get(j)+" "+targetSentence.get(i)+" "+dProbabilityMap.getCount(j, new Pair<Pair<Integer, Integer>, Integer>(new Pair<Integer, Integer>(i,I),j)));
				if(Double.isNaN(probability)) {
					System.out.println(sourceSentence.get(j)+"|"+targetSentence.get(i)+" = "+tProbabilityMap.getCount(sourceSentence.get(j), targetSentence.get(i)));
					System.out.println("d[bucket(i-j*I/J)] = "+d[bucket(i-j*I/J)]);
				}

			}
		}
		for(int i=0;i<I;i++) {
			if(alignment.getAlignedSources(i).size() > 0) {
				probability += Math.log(nProbabilityMap.getCount(alignment.getAlignedSources(i).size(),targetSentence.get(i))) + Factorial.get(alignment.getAlignedSources(i).size());
				if(Double.isNaN(probability)) {
					System.out.println(alignment.getAlignedSources(i).size()+"|"+targetSentence.get(i)+" = "+nProbabilityMap.getCount(alignment.getAlignedSources(i).size(),targetSentence.get(i)));
				}
			}
			else phi0++;
		}
		//		System.out.println("French "+sourceSentence+" English "+targetSentence);
		//		System.out.println(J+" "+phi0);
		if(Double.isNaN(probability)) {
			System.out.println("Trouble!");
		}
		if(J-2*phi0>=0) {
			double comb = 0.0;
			for(int i=1;i<=phi0;i++)
				comb += Math.log(J-phi0-i+1) - Math.log(i);
			probability += comb;
			//			probability += Math.log(Factorial.get(J - phi0));
			//			if(Double.isNaN(probability)) {
			//				System.out.println("1. J,phi0 "+J+" "+phi0);
			//				return probability;
			//			}		
			//			probability -= (Math.log(Factorial.get(J - 2*phi0)) + Math.log(Factorial.get(phi0)));
			//			if(Double.isNaN(probability)) {
			//				System.out.println("2. J,phi0 "+Factorial.get(J - 2*phi0)+" "+Factorial.get(phi0));
			//				return probability;
			//			}		
			//			probability += (J - 2*phi0)*Math.log(1 - p) + (phi0)*Math.log(p);
			if(Double.isNaN(probability)) {
				System.out.println("3. J,phi0 "+J+" "+phi0+" "+comb);
				return probability;
			}		
		}
		else return 0;
		if(Double.isNaN(probability)) {
			System.out.println("Trouble! 2");
		}
		if(Double.isNaN(Math.exp(probability))) {
			System.out.println("Trouble! exp");
		}

		return Math.exp(probability);
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#getProbSourceGivenTarget()
	 */
	@Override
	public CounterMap<String, String> getProbSourceGivenTarget() {
		return tProbabilityMap;
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#train(java.util.List)
	 */
	@Override
	public void train(List<SentencePair> trainingPairs) {

		Model2WordAligner model2 = new Model2WordAligner();
		model2.train(trainingPairs);
		CounterMap<String, String> model2ProbabilityMap = model2.getProbSourceGivenTarget();
		nCountMap = new CounterMap<String, Integer>();
		tCountMap = new CounterMap<String, String>();

		nProbabilityMap = new CounterMap<Integer, String>();
		tProbabilityMap = new CounterMap<String, String>();

		d = new double[nbuckets];
		dupdate = new double[nbuckets];

		alignmentProbabilities = new CounterMap<List<String>, List<String>>();

		phimax = 4;
		double betamax = 40000;
		Partitions = new HashMap<Integer, List<List<Integer>>>();
		Factorial = new HashMap<Integer, Integer>();
		CalculatePartitions();
		//		System.out.println(Partitions.get(4)+" "+Factorial.get(4));
		// Initialization
		p = 0.5;
		int nalignments = 100;
		possibleAlignments = new HashMap<Pair<List<String>,List<String>>, Counter<Alignment>>();
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			int I = targetWords.size();
			int J = sourceWords.size();
			CounterMap<Integer, Integer> alignmentProbability = model2.getProbAGivenFE(pair); // french, english
			Counter<Alignment> alignments = model2.CalculateTopKAlignments(alignmentProbability, nalignments, pair);
			possibleAlignments.put(new Pair<List<String>, List<String>>(sourceWords, targetWords), alignments);
			//			System.out.println("Calculated Possible Alignments");
			//			System.out.println("French "+sourceWords+" English "+targetWords);
			//			System.out.println(alignments);
			for(int i=0;i<I;i++) {
				for(int j=0;j<J;j++) {
					tProbabilityMap.setCount(sourceWords.get(j), targetWords.get(i), model2ProbabilityMap.getCount(sourceWords.get(j), targetWords.get(i)));
					//					if(alignmentProbability.getCount(j, i)>0)
					//						System.out.println("****FLAG***\n"+i+","+I+","+J+","+j);
					d[bucket(i-j*I/J)] += alignmentProbability.getCount(j, i);
				}
//				for(int phi = 1; phi <= phimax; phi++) {
//					nCountMap.setCount(targetWords.get(i), phi, 1.0);
//				}
			}

			CounterMap<Integer, Integer> alpha = new CounterMap<Integer, Integer>();
			for(int i=0;i<I;i++) {
				for(int k=1;k<=phimax;k++) {
					double beta = 0.0;
					double tempprob;
					for(int j=0;j<J;j++) {
						tempprob = model2ProbabilityMap.getCount(sourceWords.get(j), targetWords.get(i));
						//						if(tempprob == 1) {
						//							beta = betamax;
						//							break;
						////							System.out.println("Bad");
						//						}
						beta += Math.pow(tempprob/(1 - tempprob), k);
					}
					alpha.setCount(i, k, Math.pow(-1, k+1)*beta/k);
				}
			}
			for(int i=0;i<I;i++) {
				double r = 1.0;
				for(int j=0;j<J;j++) {
					r *= 1 - model2ProbabilityMap.getCount(sourceWords.get(j), targetWords.get(i));
				}
				for(int phi = 1; phi <= phimax; phi++) {
					double sum = 0.0;
					for(List<Integer> p : Partitions.get(phi)) {
						double prod = 1.0;
						for(int k=1;k<=phi;k++) {
							int n = times(p,k);
							if(n>0) {
								prod *= Math.pow(alpha.getCount(i, k),n)/Factorial.get(n);
								//								if(Double.isInfinite(alpha.getCount(i, k)))
								//									System.out.println("****FLAG***\n"+targetWords.get(i)+","+n);
							}
						}
						sum += prod;
						//						if(targetWords.get(i).equals("A")) 
						//							System.out.println(targetWords.get(i)+" += "+r+"*"+sum);
					}
					//					System.out.println(targetWords.get(i)+" += "+r+"*"+sum);
					nCountMap.incrementCount(targetWords.get(i), phi, r*sum);
				}
			}
		}
		double total;
		double dsum = 0.0;
		for(int j=0;j<nbuckets;j++)
			dsum+=d[j];
		for(int j=0;j<nbuckets;j++)
			d[j]/=dsum;

		for(String e : nCountMap.keySet()) {
			total = nCountMap.getCounter(e).totalCount();
			System.out.println(e+" "+total);

			for(Integer j : nCountMap.getCounter(e).keySet()) {
				if(total!=0)
					nProbabilityMap.setCount(j, e, nCountMap.getCount(e, j)/total);
				else
					nProbabilityMap.setCount(j, e, 0);
				nCountMap.setCount(e, j, 0);
			}
		}			
		//		for(int j=0;j<nbuckets;j++) {
		//			System.out.println(j+": "+d[j]);
		//		}
		//
		//		System.out.println(nProbabilityMap);
		//		System.out.println(tProbabilityMap);
		//		System.out.println(dProbabilityMap);

		//EM
		System.out.println("Starting EM");
		int niterations = 300;
		boolean terminate = false;
		double tolerance = 0.01;
		int q=0;
		boolean tobreak = false;

		for(q=0;q<niterations;q++) {
			System.out.println("Iteration "+q);
			if(terminate)
				break;
			terminate = true;

			//E Step
			double numerator = 0.0;
			double denominator = 0.0;
			for(int j=0;j<nbuckets;j++) {
				dupdate[j] = 0;
			}
			for(SentencePair pair : trainingPairs){
				List<String> targetWords = pair.getEnglishWords();
				List<String> sourceWords = pair.getFrenchWords();
				int I = targetWords.size();
				int J = sourceWords.size();
				Counter<Alignment> currentAlignments = possibleAlignments.get(new Pair<List<String>, List<String>>(sourceWords, targetWords));
				double probability;
				for(Alignment alignment : currentAlignments.keySet()) {
					probability = getAlignmentProbHelper(targetWords, sourceWords, alignment);
					if(Double.isNaN(probability)) {
						System.out.println("*****FLAG******\nFrench "+sourceWords+"\nEnglish: "+targetWords+"\n"+alignment);
						terminate = true;
						tobreak = true;
						break;
					}
					if(tobreak)
						break;
					//					System.out.println("Probability : "+probability);
					int phi0 = 0;
					for(int j=0;j<J;j++) {
						int i = alignment.getAlignedTarget(j);
						if(i>=0) {
							tCountMap.incrementCount(targetWords.get(i), sourceWords.get(j),probability);
							dupdate[bucket(i-j*I/J)] +=  probability;
						}
						else phi0++;
					}
					for(int i=0;i<I;i++) {
						//						if(targetWords.get(i).equals("F"))
						//							System.out.println("French "+sourceWords+" English: "+targetWords+"\n"+alignment+" "+alignment.getAlignedSources(i).size());
						if(alignment.getAlignedSources(i).size() > 0)
							nCountMap.incrementCount(targetWords.get(i), alignment.getAlignedSources(i).size(), probability);
					}
					if(J-2*phi0 >=0) {
						numerator+=phi0*probability;
						denominator+=(J-phi0)*probability;
					}
				}
				if(tobreak)
					break;

			}
			if(tobreak)
				break;
			//M Step
			p = numerator/denominator;
			dsum = 0.0;
			for(int j=0;j<nbuckets;j++)
				dsum+=dupdate[j];
			for(int j=0;j<nbuckets;j++)
				d[j]=dupdate[j]/dsum;
			//			double total;
			for(String e : nCountMap.keySet()) {
				total = nCountMap.getCounter(e).totalCount();
				System.out.println(e+" "+total);

				for(Integer j : nCountMap.getCounter(e).keySet()) {
					double oldval = nProbabilityMap.getCount(j, e);
					double num = nCountMap.getCount(e, j);
					if(num > 0) {
						nProbabilityMap.setCount(j, e, num/total);
						if(terminate && Math.abs(oldval - num/total) > tolerance) {
							terminate = false;
						}					
					}
					else {
						nProbabilityMap.setCount(j, e, 0);
						if(terminate && Math.abs(oldval) > tolerance) {
							terminate = false;
						}
					}

					nCountMap.setCount(e, j, 0);
				}
			}			
			for(String target : tCountMap.keySet()) {
				total = tCountMap.getCounter(target).totalCount();
				for(String source : tCountMap.getCounter(target).keySet()) {
					double num = tCountMap.getCount(target, source);
					double oldval = tProbabilityMap.getCount(source, target);
					if(num == 0.0) {
						if(terminate && Math.abs(oldval) > tolerance) {
							terminate = false;
						}
						tProbabilityMap.setCount(source, target, 0);
					}
					else {
						if(terminate && Math.abs(oldval - num/total) > tolerance) {
							terminate = false;
						}
						tProbabilityMap.setCount(source, target, num/total);
					}
					tCountMap.setCount(target, source, 0.0);
				}
			}

		}
		System.out.println("Converged after "+q+" iterations");
		if(tobreak)
			return;
		//		for(int j=0;j<nbuckets;j++) {
		//			System.out.println(j+": "+d[j]);
		//		}
		//
		//		System.out.println(nProbabilityMap);
		//		System.out.println(tProbabilityMap);
		//		System.out.println(dProbabilityMap);
		//		System.out.println("Calculating Alignment Probabilities");
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			//			System.out.println("**French "+sourceWords+" English "+targetWords);
			Counter<Alignment> currentAlignments = possibleAlignments.get(new Pair<List<String>, List<String>>(sourceWords, targetWords));
			total = 0.0;
			for(Alignment alignment : currentAlignments.keySet()) {
				total += getAlignmentProbHelper(targetWords, sourceWords, alignment);
				//				System.out.println(alignment+" "+total);
			}
			//			System.out.println("****French "+sourceWords+" English "+targetWords+" "+total);
			alignmentProbabilities.setCount(sourceWords, targetWords, total);
		}

	}

	void CalculatePartitions() {
		Factorial.put(0, 1);
		List<Integer> temp = new LinkedList<Integer>();
		temp.add(1);
		List<List<Integer>> temp2 = new LinkedList<List<Integer>>();
		temp2.add(temp);
		Partitions.put(1, temp2);
		Factorial.put(1, 1);
		for(int i=2;i<=phimax;i++) {
			List<List<Integer>> toadd = new LinkedList<List<Integer>>();
			List<Integer> toadd2 = new LinkedList<Integer>();
			toadd2.add(i);
			toadd.add(toadd2);
			for(int j=1;j<i;j++) {
				temp2 = Partitions.get(i-j);
				for(List<Integer> innerlist : temp2) {
					List<Integer> toadd3 = new LinkedList<Integer>();
					toadd3.add(j);
					toadd3.addAll(innerlist);
					toadd.add(toadd3);
				}
			}
			Partitions.put(i, toadd);
			Factorial.put(i, i*Factorial.get(i-1));
		}
		for(int i=phimax+1;i<Math.pow(phimax, 4);i++) {
			Factorial.put(i, i*Factorial.get(i-1));			
		}
	}

	int times(List<Integer> partition, int k) {
		int count = 0;
		for(Integer n : partition) {
			if(n==k)
				count++;
		}
		return count;
	}

	private int bucket(double x) {
		x = Math.abs(x);
		x = x/4 ; 
		if(x>= nbuckets-1)
			return nbuckets-1;
		return (int) Math.floor(x);
	}


}
