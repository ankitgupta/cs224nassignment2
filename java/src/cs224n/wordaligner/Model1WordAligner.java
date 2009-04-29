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
public class Model1WordAligner extends WordAligner {

	CounterMap<String, String> CountMap; //Target, Source
	CounterMap<String, String> ProbabilityMap; //Source, Target
	CounterMap<String, String> TempProbabilityMap; //Targer, Source
	private Counter<String> FrenchCount, EnglishCount;
//	private double epsilon = 0.1; //Probability of choosing french sentence of a particular length 
	private double NULL_PROB = 0.3; //Probability of aligning to NULL 
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
			for (int englishPosition = numEnglishWords-1; englishPosition >= 0; englishPosition--) {
				value = ((1.0 - NULL_PROB)/numEnglishWords)*ProbabilityMap.getCount(sentencePair.getFrenchWords().get(frenchPosition), sentencePair.getEnglishWords().get(englishPosition));
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
				sum+=((1.0 - NULL_PROB)/numEnglishWords)*ProbabilityMap.getCount(sourceSentence.get(frenchPosition), targetSentence.get(englishPosition));
			}
			sum+=NULL_PROB*ProbabilityMap.getCount(sourceSentence.get(frenchPosition), NULL_WORD);
			product+=Math.log(sum);
		}
		
		double probability = 0;
		for(Pair<Integer, Integer> p : alignment.getSureAlignments()) {
			if(p.getFirst()>=0) 
				probability += Math.log(((1.0 - NULL_PROB)/numEnglishWords)*ProbabilityMap.getCount(sourceSentence.get(p.getSecond()), targetSentence.get(p.getFirst())));
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

		
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			for(String source : sourceWords){
				FrenchCount.incrementCount(source, 1.0);
				for(String target : targetWords){
					CountMap.setCount(target, source, Math.random());
				}
			}
			for(String target : targetWords){
				EnglishCount.incrementCount(target, 1.0);	            
			}
		}
		double total;
		
		// Uniform
		for(String target : CountMap.keySet()) {
			total = CountMap.getCounter(target).totalCount();
			for(String source : CountMap.getCounter(target).keySet()) {
				ProbabilityMap.setCount(source, target, CountMap.getCount(target, source)/total);
				CountMap.setCount(target, source, 0.0);
			}
		}
		total = FrenchCount.keySet().size();
		for(String source : FrenchCount.keySet()) {
			ProbabilityMap.setCount(source, NULL_WORD, 1.0/total);
			CountMap.setCount(NULL_WORD, source, 0.0);
		}
		
		// Stupid
//		for(String target : CountMap.keySet()) {
//			for(String source : CountMap.getCounter(target).keySet()) {
//				TempProbabilityMap.incrementCount(target, source, Math.pow(CountMap.getCount(target, source)/(FrenchCount.getCount(source)*EnglishCount.getCount(target)),2.0));
//			}
//			total = TempProbabilityMap.getCounter(target).totalCount();
//			for(String source : CountMap.getCounter(target).keySet()) {
//				ProbabilityMap.setCount(source, target, TempProbabilityMap.getCount(target, source)/total);
//				CountMap.setCount(target, source, 0.0);				
//			}			
//		}
//		total = FrenchCount.keySet().size();
//		for(String source : FrenchCount.keySet()) {
//			ProbabilityMap.setCount(source, NULL_WORD, 1.0/total);
//			CountMap.setCount(NULL_WORD, source, 0.0);
//		}

		
		int niterations = 0;
		boolean terminate = false;
		double tolerance = 0.0005;
		int i;
		for(i=0;i<niterations;i++) {
//			for(String target : CountMap.keySet()) {
//				for(String source : CountMap.getCounter(target).keySet()) {
//					System.out.println(source+" "+target+":"+ProbabilityMap.getCount(source, target));
//				}
//			}
			if(terminate)
				break;
			terminate = true;
			// E Step
			for(SentencePair pair : trainingPairs){
				int numFrenchWords = pair.getFrenchWords().size();
				int numEnglishWords = pair.getEnglishWords().size();
				double denominator = 0;
				for (int frenchPosition = 0; frenchPosition < numFrenchWords; frenchPosition++) {
					String source = pair.getFrenchWords().get(frenchPosition);
					denominator = 0;
					for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
						denominator+=((1.0 - NULL_PROB)/numEnglishWords)*ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition));
					}
					denominator+=NULL_PROB*ProbabilityMap.getCount(source, NULL_WORD);
					for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++) {
						CountMap.incrementCount(pair.getEnglishWords().get(englishPosition), source, (((1.0 - NULL_PROB)/numEnglishWords)*ProbabilityMap.getCount(source, pair.getEnglishWords().get(englishPosition))/denominator));
					}
					CountMap.incrementCount(NULL_WORD, source, (NULL_PROB*ProbabilityMap.getCount(source, NULL_WORD)/denominator));
				}
			}
			// M Step
			for(String target : CountMap.keySet()) {
				total = CountMap.getCounter(target).totalCount();
				for(String source : CountMap.getCounter(target).keySet()) {
					double oldval = ProbabilityMap.getCount(source, target);
					if(terminate && Math.abs(oldval - CountMap.getCount(target, source)/total) > tolerance) {
						terminate = false;
					}
					ProbabilityMap.setCount(source, target, CountMap.getCount(target, source)/total);
					CountMap.setCount(target, source, 0.0);
				}
			}
		}
		System.out.println("Terminated after "+i+" iterations");
	
		
	}

}
