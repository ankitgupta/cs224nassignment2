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
public class StupidWordAligner extends WordAligner {

	private CounterMap<String,String> JointCount;
	private CounterMap<String,String> ProbabilityMap, FinalProbabilityMap;
	private Counter<String> FrenchCount, EnglishCount;
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
				value = FinalProbabilityMap.getCount(sentencePair.getFrenchWords().get(frenchPosition), sentencePair.getEnglishWords().get(englishPosition));
				if(max < value) {
					max = value;
					index = englishPosition;
				}
			}
			value = FinalProbabilityMap.getCount(sentencePair.getFrenchWords().get(frenchPosition), NULL_WORD);
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
		double probability = 0;
		for(Pair<Integer, Integer> p : alignment.getSureAlignments()) {
			if(p.getFirst()>=0) 
				probability += Math.log(FinalProbabilityMap.getCount(sourceSentence.get(p.getSecond()), targetSentence.get(p.getFirst())));
			else
				probability += Math.log(FinalProbabilityMap.getCount(sourceSentence.get(p.getSecond()), NULL_WORD));
		}
		return Math.exp(probability);
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#getProbSourceGivenTarget()
	 */
	@Override
	public CounterMap<String, String> getProbSourceGivenTarget() {
		return FinalProbabilityMap;
	}

	/* (non-Javadoc)
	 * @see cs224n.wordaligner.WordAligner#train(java.util.List)
	 */
	@Override
	public void train(List<SentencePair> trainingPairs) {
		JointCount = new CounterMap<String, String>();
		FrenchCount = new Counter<String>();
		EnglishCount = new Counter<String>();
		ProbabilityMap = new CounterMap<String, String>();
		FinalProbabilityMap = new CounterMap<String, String>();
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			for(String source : sourceWords){
				FrenchCount.incrementCount(source, 1.0);
				for(String target : targetWords){
					JointCount.incrementCount(target, source, 1.0);
				}
			}
			for(String target : targetWords){
				EnglishCount.incrementCount(target, 1.0);	            
			}
		}
		for(String target : JointCount.keySet()) {
			for(String source : JointCount.getCounter(target).keySet()) {
				ProbabilityMap.incrementCount(target, source, Math.pow(JointCount.getCount(target, source),3.0)/(FrenchCount.getCount(source)*EnglishCount.getCount(target)));
//				ProbabilityMap.incrementCount(target, source, JointCount.getCount(target, source)/(EnglishCount.getCount(target)));
//				ProbabilityMap.incrementCount(target, source, Math.pow(JointCount.getCount(target, source),2));
			}
			double total = ProbabilityMap.getCounter(target).totalCount();
			for(String source : JointCount.getCounter(target).keySet()) {
				FinalProbabilityMap.setCount(source, target, ProbabilityMap.getCount(target, source)/total);
//				System.out.println(target+" "+source+":"+FinalProbabilityMap.getCount(source, target));
			}			
		}
		double total = FrenchCount.keySet().size();
		for(String source : FrenchCount.keySet()) {
			FinalProbabilityMap.setCount(source, NULL_WORD, 1.0/total);
//			System.out.println(NULL_WORD+" "+source+":"+FinalProbabilityMap.getCount(source, NULL_WORD));			
		}
	}

}
