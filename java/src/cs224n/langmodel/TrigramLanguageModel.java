/**
 * 
 */
package cs224n.langmodel;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.omg.CORBA.UNKNOWN;

import cs224n.util.Counter;
import cs224n.util.CounterMap;

/**
 * @author ankit
 *
 */
public class TrigramLanguageModel implements LanguageModel {
	
	public static final String STOP = "</S>";
	public static final String UNKNOWN = "**UNKNOWN**";
	
	HashMap<String, CounterMap<String, String>> trigramCounter;
	CounterMap<String, String> bigramCounter;
	Counter<String> unigramCounter;
	CounterMap<String, String> ADNZBigramCounter;
	CounterMap<String, String> normalizationCounter;
	Counter<String> bigramnormalizationCounter;
	
	CounterMap<String, String> KNBigramCounter;
	Counter<String> KNUnigramCounter;
	double KNUnigramtotal;
	
	CounterMap<String, String> NZBigramCounter;
	Counter<String> NZUnigramCounter;
	
	double total;
	private Counter<String> sumPwi;
	
	double lambda1=0.1, lambda2=0.2, ideallambda1, ideallambda2, min=10000000.0;
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#checkModel()
	 */
	public double checkModel() {
		// TODO Auto-generated method stub
		int nsamples = 50;
		double tolerance = 0.005;
		String firstword;
		Iterator<String> iter = trigramCounter.keySet().iterator();
		for(int i=0; i<nsamples; i++) {
			firstword = iter.next();
			double probability = 0.0;
			String secondword = trigramCounter.get(firstword).keySet().iterator().next();
			
			for(String word : unigramCounter.keySet()) {
				probability += getWordProbability(firstword, secondword, word);
			}
			probability += getWordProbability(firstword, secondword, UNKNOWN);
			if(Math.abs(1 - probability) > tolerance)
				return probability;
		}
		return 1.0;
	}
	
	/**
	 * Returns a random sentence sampled according to the model.  We generate
	 * words until the stop token is generated, and return the concatenation.
	 */
	public List<String> generateSentence() {
		List<String> sentence = new ArrayList<String>();
		String firstword = STOP, secondword = STOP;
		
		String word = generateWord(firstword, secondword);
		while (!word.equals(STOP)) {
			sentence.add(word);
			firstword = secondword;
			secondword = word;
			word = generateWord(firstword,secondword);
		}
		return sentence;
	}
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#generateSentence()
	 */
	private String generateWord(String firstword, String secondword) {
		double sample = Math.random();
		double sum = 0.0;
		for (String word : unigramCounter.keySet()) {
			sum += getWordProbability(firstword, secondword, word);
			if (sum > sample) {
				//				System.out.println(firstword+","+secondword+","+word+":"+getWordProbability(firstword, secondword, word));
				return word;
			}
		}
		return UNKNOWN;   // a little probability mass was reserved for unknowns
	}
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#getSentenceProbability(java.util.List)
	 */
	public double getSentenceProbability(List<String> sentence) {
		List<String> stoppedSentence = new ArrayList<String>(sentence);
		stoppedSentence.add(STOP);
		
		double probability = 0.0;
		for (int index = 0; index < stoppedSentence.size(); index++) {
			probability += Math.log(getWordProbability(stoppedSentence, index));
		}
		return Math.exp(probability);
	}
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#getUnigramCounter()
	 */
	public Counter<String> getUnigramCounter() {
		// TODO Auto-generated method stub
		return null;
	}
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#getWordCounter()
	 */
	public CounterMap<String, String> getBigramCounter() {
		// TODO Auto-generated method stub
		return bigramCounter;
	}
	
	public 	HashMap<String, CounterMap<String, String>> getTrigramCounter() {
		return trigramCounter;
	}
	
	//Bigram Probability
	private double getWordProbability(String w1, String w2) {
		//		return getWordProbabilityCustomSmoothing(w1, w2);
		//		return getWordProbabilityLaplaceSmoothing(w1, w2);
		return getWordProbabilityAbsoluteDiscounting(w1, w2);
		//		return getWordProbabilityKneserNey(w1, w2);
		//		return getWordProbabilityKneserNeyBackoff(w1, w2);
		
	}
	public double getWordProbabilityCustomSmoothing(String w1, String w2) {
		
		// CustomSmoothing
		
		double bigramcount = bigramCounter.getCount(w1, w2);
		double probability = 0.0;
		
		if (bigramCounter.keySet().contains(w1)) {
			double numBigramsWithZeroCount = unigramCounter.size()
			- bigramCounter.getCounter(w1).size();
			if (bigramcount > 0) {
				probability = ((double) bigramcount)
				/ (1.0 + (double) bigramnormalizationCounter.getCount(w1));
			} else {
				probability = (1.0 / (1 + numBigramsWithZeroCount))
				/ (1.0 + (double) bigramnormalizationCounter.getCount(w1));
			}
			return probability;
		} else {
			if ((unigramCounter.getCount(w2) == 0)) {
				return 1 / (total + 1.0);
			}
			return unigramCounter.getCount(w2) / (1.0 + total);
		}
	}
	
	public double getWordProbabilityLaplaceSmoothing(String w1, String w2) {
		
		// LaplaceSmoothing
		
		double probability = 0.0;
		
		if (bigramCounter.keySet().contains(w1)) {
			double bigramcount = bigramCounter.getCount(w1, w2);
			double vocab = unigramCounter.size();
			double numBigramsWithZeroCount = unigramCounter.size()
			- bigramCounter.getCounter(w1).size();
			if (bigramcount > 0) {
				probability = ((double) bigramcount + 1)
				/ (1.0 + vocab + (double) bigramnormalizationCounter
				   .getCount(w1));
			} else {
				probability = (1.0) / (1.0 + vocab + (double) bigramnormalizationCounter
									   .getCount(w1));
			}
			return probability;
		} else {
			if ((unigramCounter.getCount(w2) == 0)) {
				return 1 / (total + 1.0);
			}
			return unigramCounter.getCount(w2) / (1.0 + total);
			
		}
	}
	
	public double getWordProbabilityAbsoluteDiscounting(String w1, String w2) {
		
		// AbsoluteDiscounting
		
		double probability = 0.0;
		double discount = 0.57;
		if (bigramCounter.keySet().contains(w1)) {
			double bigramcount = bigramCounter.getCount(w1, w2);
			double numBigramsWithZeroCount = unigramCounter.size()
			- bigramCounter.getCounter(w1).size();
			if (bigramcount > 0) {
				probability = ((double) bigramcount - discount)
				/ ((double) bigramnormalizationCounter.getCount(w1));
			} else {
				double alpha = (discount * ((bigramCounter.getCounter(w1))
											.size()))
				/ (bigramnormalizationCounter.getCount(w1));
				alpha = alpha / sumPwi.getCount(w1);
				if (unigramCounter.getCount(w2) == 0) {
					probability = alpha * (1 / (1 + total));
				} else {
					probability = alpha
					* (unigramCounter.getCount(w2) / (1 + total));
				}
			}
			return probability;
		} else {
			if ((unigramCounter.getCount(w2) == 0)) {
				return 1 / (total + 1.0);
			}
			return unigramCounter.getCount(w2) / (1.0 + total);
		}
	}
	
	//Unigram Probability 	
	private double getWordProbability(String word) {
		//		return getWordProbabilityCustomSmoothing(word);
		//		return getWordProbabilityLaplaceSmoothing(word);
		return getWordProbabilityAbsoluteDiscounting(word);
		//		return getWordProbabilityKneserNey(word);
	}
	
	private double getWordProbabilityCustomSmoothing(String word) {
		double count = unigramCounter.getCount(word);
		
		if (count == 0) {
			return 1.0 / (total + 1.0);
		}
		return count / (total + 1.0);		
	}
	
	private double getWordProbabilityLaplaceSmoothing(String word) {
		double count = unigramCounter.getCount(word);
		
		if (count == 0) {
			return 1.0 / (total + unigramCounter.keySet().size() + 1.0);
		}
		return (count + 1.0) / (total + unigramCounter.keySet().size() + 1.0);		
	}
	private double getWordProbabilityAbsoluteDiscounting(String word) {
		double count = unigramCounter.getCount(word);
		double discount = 0.1;
		if (count == 0) {
			return discount*unigramCounter.keySet().size() / total;
		}
		return (count - discount) / total;		
	}
	
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#getWordProbability(java.util.List, int)
	 */
	public double getWordProbability(List<String> sentence, int index) {
		String firstword = STOP, secondword = STOP;
		if(index >= 1)
			secondword = sentence.get(index-1);
		if(index >= 2)
			firstword = sentence.get(index-2);
		String currentword = sentence.get(index);
		return getWordProbability(firstword, secondword, currentword);
	}
	private double getWordProbability(String firstword, String secondword, String thirdword) {
		//		return getWordProbabilityKneserNeySmoothingBackoff(firstword, secondword, thirdword);
		//		return getWordProbabilityKneserNeySmoothing(firstword, secondword, thirdword);
		//		return getWordProbabilityLaplaceSmoothing(firstword, secondword, thirdword);
		//		return getWordProbabilityDiscountSmoothing(firstword, secondword, thirdword);
		//		return getWordProbabilityCustomSmoothing(firstword, secondword, thirdword);
		return getWordProbabilityInterpolation(firstword, secondword, thirdword);
	}
	private double getWordProbabilityInterpolation(String firstword, String secondword, String thirdword) {
		return lambda1*getWordProbability(thirdword)+lambda2*getWordProbability(secondword, thirdword)+(1-lambda1-lambda2)*getWordProbabilityDiscountSmoothing(firstword, secondword, thirdword);
		
	}
	public double getWordProbabilityLaplaceSmoothing(String firstword, String secondword, String thirdword) {
		// Laplace Smoothing
		double numerator;
		if(!trigramCounter.containsKey(firstword)) {
			return getWordProbability(secondword, thirdword);
		}
		if(!trigramCounter.get(firstword).keySet().contains(secondword)) {
			if(!unigramCounter.containsKey(secondword))
				return getWordProbability(thirdword);
			else
				return getWordProbability(secondword, thirdword);
		}
		numerator = trigramCounter.get(firstword).getCount(secondword, thirdword);
		double denominator = normalizationCounter.getCount(firstword, secondword);
		if(numerator == 0)
			return 1.0 / (denominator + unigramCounter.keySet().size() + 1.0);
		
		return (numerator + 1.0) / (denominator + unigramCounter.keySet().size() + 1.0);
	}
	
	private double getWordProbabilityKneserNey(String word) {
		double count = KNUnigramCounter.getCount(word);
		if (count == 0) {
			return 1.0 / (KNUnigramtotal + 1.0);
		}
		return count / (KNUnigramtotal + 1.0);
		
	}
	
	private double getWordProbabilityKneserNey(String firstword, String secondword) {
		double discount = 0.85; 
		if (!bigramCounter.keySet().contains(firstword)) {
			return getWordProbabilityKneserNey(secondword) ; 
		}
		double numerator = bigramCounter.getCount(firstword, secondword);
		double denominator = bigramnormalizationCounter.getCount(firstword);
		if(numerator == 0)
			return discount*bigramCounter.getCounter(firstword).keySet().size()*getWordProbabilityKneserNey(secondword)/denominator;
		return (numerator - discount) / denominator + discount*bigramCounter.getCounter(firstword).keySet().size()*getWordProbabilityKneserNey(secondword)/denominator;
	}
	
	public double getWordProbabilityKneserNeySmoothing(String firstword, String secondword, String thirdword) {
		// Kneser Ney Smoothing
		
		double discount = 0.97; 
		if(!trigramCounter.containsKey(firstword)) {
			return getWordProbabilityKneserNey(secondword, thirdword);
		}
		if(!trigramCounter.get(firstword).keySet().contains(secondword)) {
			if(!unigramCounter.containsKey(secondword))
				return getWordProbability(thirdword);
			else
				return getWordProbabilityKneserNey(secondword, thirdword);
		}
		double numerator = trigramCounter.get(firstword).getCount(secondword, thirdword);
		double denominator = normalizationCounter.getCount(firstword, secondword);
		if(numerator == 0)
			return discount*trigramCounter.get(firstword).getCounter(secondword).keySet().size()*getWordProbabilityKneserNey(secondword, thirdword) / denominator;
		
		return (numerator - discount) / denominator + discount*trigramCounter.get(firstword).getCounter(secondword).keySet().size()*getWordProbabilityKneserNey(secondword, thirdword) / denominator;
	}
	
	private double getWordProbabilityKneserNeyBackoff(String firstword, String secondword) {
		double discount = 0.85;
		if (!bigramCounter.keySet().contains(firstword)) {
			return getWordProbabilityKneserNey(secondword) ; 
		}
		double numerator = bigramCounter.getCount(firstword, secondword);
		double denominator = bigramnormalizationCounter.getCount(firstword);
		if(numerator == 0)
			return discount*bigramCounter.getCounter(firstword).keySet().size()*getWordProbabilityKneserNey(secondword)/(denominator*(1 - NZUnigramCounter.getCount(firstword)));
		return (numerator - discount) / denominator;
	}
	
	public double getWordProbabilityKneserNeySmoothingBackoff(String firstword, String secondword, String thirdword) {
		// Kneser Ney Smoothing
		double discount = 0.75;
		
		double numerator;
		if(!trigramCounter.containsKey(firstword)) {
			return getWordProbabilityKneserNey(secondword, thirdword);
		}
		if(!trigramCounter.get(firstword).keySet().contains(secondword)) {
			if(!unigramCounter.containsKey(secondword))
				return getWordProbability(thirdword);
			else
				return getWordProbabilityKneserNey(secondword, thirdword);
		}
		numerator = trigramCounter.get(firstword).getCount(secondword, thirdword);
		double denominator = normalizationCounter.getCount(firstword, secondword);
		if(numerator == 0) {
			return discount*trigramCounter.get(firstword).getCounter(secondword).keySet().size()*getWordProbabilityKneserNey(secondword, thirdword) / (denominator *(1 - NZBigramCounter.getCount(firstword,secondword)));
		}
		
		return (numerator - discount) / denominator;
	}
	
	public double getWordProbabilityDiscountSmoothing(String firstword, String secondword, String thirdword) {
		// Discount Smoothing
		
		double discount = 0.70; 
		
		double numerator;
		if(!trigramCounter.containsKey(firstword)) {
			return getWordProbability(secondword, thirdword);
		}
		if(!trigramCounter.get(firstword).keySet().contains(secondword)) {
			if(!unigramCounter.containsKey(secondword))
				return getWordProbability(thirdword);
			else
				return getWordProbability(secondword, thirdword);
		}
		numerator = trigramCounter.get(firstword).getCount(secondword, thirdword);
		double denominator = normalizationCounter.getCount(firstword, secondword);
		if(numerator == 0)
			return discount*trigramCounter.get(firstword).getCounter(secondword).keySet().size()*getWordProbability(secondword, thirdword) / ((1-ADNZBigramCounter.getCount(firstword, secondword))*denominator);
		
		return (numerator - discount) / (denominator);
	}
	
	public double getWordProbabilityCustomSmoothing(String firstword, String secondword, String thirdword) {
		double numerator;
		if(!trigramCounter.containsKey(firstword)) {
			return getWordProbability(secondword, thirdword);
		}
		if(!trigramCounter.get(firstword).keySet().contains(secondword)) {
			if(!unigramCounter.containsKey(secondword))
				return getWordProbability(thirdword);
			else
				return getWordProbability(secondword, thirdword);
		}
		numerator = trigramCounter.get(firstword).getCount(secondword, thirdword);
		double denominator = normalizationCounter.getCount(firstword, secondword);
		if(numerator == 0)
			return 1.0 / ((unigramCounter.keySet().size() - trigramCounter.get(firstword).getCounter(secondword).keySet().size() + 1.0)*(denominator + 1.0));
		
		return numerator / (denominator + 1.0);
		
	}
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#train(java.util.Collection)
	 */
	public void train(Collection<List<String>> trainingSentences) {
		// TODO Auto-generated method stub
		trigramCounter = new HashMap<String, CounterMap<String, String>>();
		bigramCounter = new CounterMap<String, String>();
		unigramCounter = new Counter<String>();
		
		KNBigramCounter = new CounterMap<String, String>();
		String firstword = STOP, secondword = STOP;
		
		for (List<String> sentence : trainingSentences) {
			List<String> stoppedSentence = new ArrayList<String>(sentence);
			stoppedSentence.add(STOP);
			firstword = STOP;
			secondword = STOP;
			bigramCounter.incrementCount(firstword, secondword, 1.0);
			
			for (String word : stoppedSentence) {
				unigramCounter.incrementCount(word, 1.0);
				bigramCounter.incrementCount(secondword, word, 1.0);
				KNBigramCounter.incrementCount(word, secondword, 1.0);
				if(!trigramCounter.containsKey(firstword)) {
					trigramCounter.put(firstword, new CounterMap<String, String>());
				}
				trigramCounter.get(firstword).incrementCount(secondword, word, 1.0);
				firstword = secondword;
				secondword = word;
			}
		}
		//For Bigram
		sumPwi = new Counter<String>();
		total = unigramCounter.totalCount() ; 
		
		for (Iterator iterator1 = unigramCounter.keySet().iterator(); iterator1.hasNext();) {
			String w1 = (String) iterator1.next();
			sumPwi.setCount(w1, 0);
		}		
		
		for (Iterator iterator1 =    bigramCounter.keySet().iterator(); iterator1.hasNext();) {
			String w1 = (String) iterator1.next();
			for (Iterator iterator2 = bigramCounter.getCounter(w1).keySet().iterator(); iterator2.hasNext();) {
				String w2= (String) iterator2.next();
				sumPwi.incrementCount(w1, unigramCounter.getCount(w2)/ ( 1 + total )); 
			}
		}
		
		for (Iterator iterator1 = unigramCounter.keySet().iterator(); iterator1.hasNext();) {
			String w1 = (String) iterator1.next();
			sumPwi.setCount(w1, 1- sumPwi.getCount(w1));
		}
		//End For Bigram
		
		
		performKneserNeySmoothing();
		computeNonZeroProbabilities();
		
		//		System.out.println(getWordProbability(STOP, STOP, "i"));
		//		double max = 0, temp;
		//		String w="";
		//		for(String word : unigramCounter.keySet()) {
		//			temp = getWordProbability("amongst", "all", word);
		//			if(temp>max) {
		//				max=temp;
		//				w=word;
		//			}
		//			
		//		}
		//		temp = getWordProbability("amongst", "all", UNKNOWN);
		//		if(temp>max) {
		//			max=temp;
		//			w=UNKNOWN;
		//		}
		//		
		//			System.out.println("amongst all "+w+":"+max);
		//		System.out.println("this week we want to update you on some further management changes that have taken shape");
		//		System.out.println("on some further:"+getWordProbability("on", "some", "further"));
		//		System.out.println("some further management:"+getWordProbability("some", "further", "management"));
		//		System.out.println("further management changes:"+getWordProbability("further", "management", "changes"));
		//		System.out.println("management changes that:"+getWordProbability("management", "changes", "that"));
		//
		//		System.out.println("\n\nthis week we want to update you on some management further changes that have taken shape");
		//		System.out.println("on some management:"+getWordProbability("on", "some", "management"));
		//		System.out.println("some management further:"+getWordProbability("some", "management", "further"));
		//		System.out.println("management further changes:"+getWordProbability("management", "further", "changes"));
		//		System.out.println("further changes that:"+getWordProbability("further", "changes", "that"));
		
		
	}
	
	private void computeNonZeroProbabilities() {
		normalizationCounter = new CounterMap<String, String>();
		for(String firstword : trigramCounter.keySet()) {
			for(String secondword : trigramCounter.get(firstword).keySet()) {
				normalizationCounter.incrementCount(firstword, secondword, trigramCounter.get(firstword).getCounter(secondword).totalCount());
			}
		}
		
		
		bigramnormalizationCounter = new Counter<String>();
		for(String firstword : bigramCounter.keySet()) {
			bigramnormalizationCounter.incrementCount(firstword, bigramCounter.getCounter(firstword).totalCount());
		}
		
		NZBigramCounter = new CounterMap<String, String>();
		NZUnigramCounter = new Counter<String>();
		double probability = 0.0;
		for(String firstword : bigramCounter.keySet()) {
			probability = 0.0;
			for(String secondword : bigramCounter.getCounter(firstword).keySet()) {
				probability += getWordProbabilityKneserNey(secondword);
			}
			NZUnigramCounter.incrementCount(firstword, probability);
		}
		for(String firstword : trigramCounter.keySet()) {
			for(String secondword : trigramCounter.get(firstword).keySet()) {
				probability = 0.0;
				for(String thirdword : trigramCounter.get(firstword).getCounter(secondword).keySet()) {
					probability += getWordProbabilityKneserNey(secondword, thirdword);
				}
				NZBigramCounter.incrementCount(firstword, secondword, probability);
			}
		}
		
		ADNZBigramCounter = new CounterMap<String, String>();
		for(String firstword : trigramCounter.keySet()) {
			for(String secondword : trigramCounter.get(firstword).keySet()) {
				probability = 0.0;
				for(String thirdword : trigramCounter.get(firstword).getCounter(secondword).keySet()) {
					probability += getWordProbability(secondword, thirdword);
				}
				ADNZBigramCounter.incrementCount(firstword, secondword, probability);
			}
		}
		
		
	}
	
	private void performKneserNeySmoothing(){ 
		KNUnigramCounter = new Counter<String>();
		for(String word : unigramCounter.keySet()) {
			KNUnigramCounter.incrementCount(word, KNBigramCounter.getCounter(word).keySet().size()/KNBigramCounter.getCounter(word).totalCount());
		}
		KNUnigramtotal = KNUnigramCounter.totalCount();
		System.out.println("Done Unigram Keneser Ney Smoothing");
	}
	
	
	/* (non-Javadoc)
	 * @see cs224n.langmodel.LanguageModel#updateWordProbability(java.lang.String, double)
	 */
	public void updateWordProbability(String word, double probability) {
		// TODO Auto-generated method stub
		
	}
	
	public boolean validate(double perplexity) {
		if(perplexity<min) {
			min = perplexity;
			ideallambda1 = lambda1;
			ideallambda2 = lambda2;
		}
		double tolerance = 0.05;
	    NumberFormat nf = new DecimalFormat("0.0000");
		
		if(Math.abs(lambda1-1.0)<=tolerance) {
			lambda1 = ideallambda1;
			lambda2 = ideallambda2;
			System.out.println("*Best* Unigram weight: "+nf.format(lambda1)+" Bigram weight:"+nf.format(lambda2)+" Trigram weight:"+nf.format(1.0-lambda1-lambda2));
			
			return true;
		}
		else {
			if(Math.abs(lambda2+lambda1-1.0)<=tolerance) {
				lambda2=0.0;
				lambda1+=0.1;
			}
			else
				lambda2+=0.1;
		}
		System.out.println("Unigram weight: "+nf.format(lambda1)+" Bigram weight:"+nf.format(lambda2)+" Trigram weight:"+nf.format(1.0-lambda1-lambda2));
		return false;
	}
	
}
