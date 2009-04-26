/**
 * 
 */
package cs224n.wordaligner;

import java.util.List;

import cs224n.util.Alignment;
import cs224n.util.Counter;
import cs224n.util.CounterMap;
import cs224n.util.Pair;
import cs224n.util.SentencePair;

/**
 * @author ankit
 *
 */
public class Model2WordAligner extends WordAligner {

	CounterMap<String, String> CountMap; //Target, Source
	CounterMap<String, String> ProbabilityMap; //Source, Target
	CounterMap<String, String> TempProbabilityMap; //Targer, Source
	private Counter<String> FrenchCount, EnglishCount;
//	private double epsilon = 0.1; //Probability of choosing french sentence of a particular length 
	private double NULL_PROB = 0.2; //Probability of aligning to NULL 
	private int nbuckets = 300;
	private double[] d, dupdate;
	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#alignSentencePair(cs224n.util.SentencePair)
	 */
	@Override
	public Alignment alignSentencePair(SentencePair sentencePair) {
		// TODO Auto-generated method stub
		Alignment alignment = new Alignment();
		int numFrenchWords = sentencePair.getFrenchWords().size();
		int numEnglishWords = sentencePair.getEnglishWords().size();
		double max=0.0, value;
		int index = -1;
		for (int frenchPosition = 0; frenchPosition < numFrenchWords; frenchPosition++) {
			max = 0.0;
			for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
				value = d[bucket(englishPosition - (double)frenchPosition*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(sentencePair.getFrenchWords().get(frenchPosition), sentencePair.getEnglishWords().get(englishPosition));
				if(max < value) {
					max = value;
					index = englishPosition;
				}
			}
			value = NULL_PROB*ProbabilityMap.getCount(sentencePair.getFrenchWords().get(frenchPosition), NULL_WORD);
			if(max < value) {
				max = value;
				index = -1;
			}
			alignment.addAlignment(index, frenchPosition, true);
		}
		return alignment;
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#getAlignmentProb(java.util.List, java.util.List, cs224n.util.Alignment)
	 */
	@Override
	public double getAlignmentProb(List<String> targetSentence,
			List<String> sourceSentence, Alignment alignment) {
		// TODO Auto-generated method stub
		
		int numFrenchWords = sourceSentence.size();
		int numEnglishWords = targetSentence.size();
		double sum=0;
		double product=0;
		for (int frenchPosition = 0; frenchPosition < numFrenchWords; frenchPosition++) {
			sum=0;
			for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
				sum+=d[bucket(englishPosition - (double)frenchPosition*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(sourceSentence.get(frenchPosition), targetSentence.get(englishPosition));
			}
			sum+=NULL_PROB*ProbabilityMap.getCount(sourceSentence.get(frenchPosition), NULL_WORD);
			product+=Math.log(sum);
		}
		
		double probability = 0;
		for(Pair<Integer, Integer> p : alignment.getSureAlignments()) {
			if(p.getFirst()>=0) 
				probability += Math.log(d[bucket(p.getFirst() - (double)p.getSecond()*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(sourceSentence.get(p.getSecond()), targetSentence.get(p.getFirst())));
			else
				probability += Math.log(NULL_PROB*ProbabilityMap.getCount(sourceSentence.get(p.getSecond()), NULL_WORD));
		}
		return Math.exp(probability-product);
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#getProbSourceGivenTarget()
	 */
	@Override
	public CounterMap<String, String> getProbSourceGivenTarget() {
		// TODO Auto-generated method stub
		return ProbabilityMap;
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#train(java.util.List)
	 */
	@Override
	public void train(List<SentencePair> trainingPairs) {
		// TODO Auto-generated method stub
		CountMap = new CounterMap<String, String>();
		ProbabilityMap = new CounterMap<String, String>();
		TempProbabilityMap = new CounterMap<String, String>();
		
		FrenchCount = new Counter<String>();
		EnglishCount = new Counter<String>();

		d = new double[nbuckets];
		dupdate = new double[nbuckets];
		
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			for(String source : sourceWords){
				FrenchCount.incrementCount(source, 1.0);
				for(String target : targetWords){
					CountMap.setCount(target, source, 1.0);
				}
			}
			for(String target : targetWords){
				EnglishCount.incrementCount(target, 1.0);	            
			}
		}
		double total;
		
		// Uniform
//		for(String target : CountMap.keySet()) {
//			total = CountMap.getCounter(target).totalCount();
//			for(String source : CountMap.getCounter(target).keySet()) {
//				ProbabilityMap.setCount(source, target, 1.0/total);
//				CountMap.setCount(target, source, 0.0);
//			}
//		}
//		total = FrenchCount.keySet().size();
//		for(String source : FrenchCount.keySet()) {
//			ProbabilityMap.setCount(source, NULL_WORD, 1.0/total);
//			CountMap.setCount(NULL_WORD, source, 0.0);
//		}
		
		// Stupid
		for(String target : CountMap.keySet()) {
			for(String source : CountMap.getCounter(target).keySet()) {
				TempProbabilityMap.incrementCount(target, source, Math.pow(CountMap.getCount(target, source)/(FrenchCount.getCount(source)*EnglishCount.getCount(target)),2.0));
			}
			total = TempProbabilityMap.getCounter(target).totalCount();
			for(String source : CountMap.getCounter(target).keySet()) {
				ProbabilityMap.setCount(source, target, TempProbabilityMap.getCount(target, source)/total);
				CountMap.setCount(target, source, 0.0);				
			}			
		}
		total = FrenchCount.keySet().size();
		for(String source : FrenchCount.keySet()) {
			ProbabilityMap.setCount(source, NULL_WORD, 1.0/total);
			CountMap.setCount(NULL_WORD, source, 0.0);
		}

		for(int i=0;i<nbuckets;i++)
			d[i] = (1.0 - NULL_PROB)/(nbuckets);
		
		int niterations = 120;
		for(int i=0;i<niterations;i++) {
//			for(String target : CountMap.keySet()) {
//				for(String source : CountMap.getCounter(target).keySet()) {
//					System.out.println(source+" "+target+":"+ProbabilityMap.getCount(source, target));
//				}
//			}

			// E Step
			for(int j=0;j<nbuckets;j++) {
//				System.out.println(i+":"+j+": "+dupdate[j]);
				dupdate[j] = 0;
			}
			for(SentencePair pair : trainingPairs){
				int numFrenchWords = pair.getFrenchWords().size();
				int numEnglishWords = pair.getEnglishWords().size();
				double denominator = 0;
				for (int frenchPosition = 0; frenchPosition < numFrenchWords; frenchPosition++) {
					String source = pair.getFrenchWords().get(frenchPosition);
					denominator = 0;
					for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
						denominator+=d[bucket(englishPosition - (double)frenchPosition*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition));
//						System.out.println(frenchPosition+","+englishPosition+":"+denominator+" "+d[bucket(englishPosition - frenchPosition*numEnglishWords/numFrenchWords)]+" "+ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition)));
					}
					denominator+=NULL_PROB*ProbabilityMap.getCount(source, NULL_WORD);
//					System.out.println(denominator);
					for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
						CountMap.incrementCount(pair.getEnglishWords().get(englishPosition), source, d[bucket(englishPosition - (double)frenchPosition*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition))/denominator);
						dupdate[bucket(englishPosition - (double)frenchPosition*numEnglishWords/numFrenchWords)] += d[bucket(englishPosition - frenchPosition*numEnglishWords/numFrenchWords)]*ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition))/denominator;
					}
					CountMap.incrementCount(NULL_WORD, source, (NULL_PROB*ProbabilityMap.getCount(source, NULL_WORD)/denominator));
				}
			}
			double dsum = 0.0;
			for(int j=0;j<nbuckets;j++)
				dsum+=dupdate[j];

			// M Step			
			for(int j=0;j<nbuckets;j++) {
				dupdate[j] = (1-NULL_PROB)*dupdate[j]/(dsum);
				d[j] = dupdate[j];
			}
			
			for(String target : CountMap.keySet()) {
				total = CountMap.getCounter(target).totalCount();
				for(String source : CountMap.getCounter(target).keySet()) {
					double numerator = CountMap.getCount(target, source);
					if(numerator == 0.0)
						ProbabilityMap.setCount(source, target, 0);
					else
						ProbabilityMap.setCount(source, target, numerator/total);
					CountMap.setCount(target, source, 0.0);
				}
			}
		}
	
		for(int j=0;j<nbuckets;j++) {
			System.out.println(j+": "+d[j]);
		}

		
	}

	private int bucket(double x) {
		x = Math.abs(x);
		if(x*4.0 >= nbuckets)
			return nbuckets-1;
		return (int) Math.floor(x*4.0);
	}
}