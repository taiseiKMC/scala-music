import java.io.File;
import java.io.StringWriter;
import java.io.PrintWriter;
import java.net.URL;
import java.util.Arrays;
import java.util.ResourceBundle;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.IntStream;

import org.apache.commons.math3.complex.Complex;

object PitchDetection {

    //自己相関
    def autocorelation(ary:Array[Double], minRate:Double, maxRate:Double, sampleRate:Double):Double =
    {
        var ma= -1000000000.0
        var maind=0
        val mint=(sampleRate/minRate).asInstanceOf[Int]
        val maxt=(sampleRate/maxRate).asInstanceOf[Int]
        for (i <- mint to maxt by -1)
        {
            var sum=0.0
            for(j <- 0 until ary.length-i)
            {
                sum+=ary(j)*ary(j+i)
            }
            if(ma<sum)
            {
                ma=sum
                maind=i
            }
        }
        if (maind==maxt) return 0.0
        return sampleRate/maind
    }

    //SHS
    //ary...spectrum
    def subharmonicSummation(ary:Array[Double], minRate:Double, maxRate:Double, sampleRate:Double):Double=
    {
        val mint=(ary.length*minRate/sampleRate).asInstanceOf[Int]
        val maxt=(ary.length*maxRate/sampleRate).asInstanceOf[Int]
        var sum=Arrays.copyOf(ary, ary.length)
        var h=0.84
        val hf=h
        var ma= -1000.0
        var maind= -1
        for (i <- mint until maxt)
        {
            h=hf
            var j=2
            while(i*j < ary.length)
            {
                sum(i)+=h*ary(i*j)
                h*=hf
                j+=1
            }

            if(ma<sum(i))
            {
                ma=sum(i)
                maind=i
            }
        }
        return maind*sampleRate/ary.length
    }

    private def nsdf(ary:Array[Double]):Array[Double] =
    {
        val r=new Array[Double](ary.length)
        for(i <- 0 until ary.length; j <- 0 until ary.length-i)
        {
            r(i)+=ary(j)*ary(j+i)
        }
        val ssum=new Array[Double](ary.length+1)
        for(i <- 0 until ary.length)
            ssum(i+1)=ssum(i) + Math.pow(ary(i),2)
        for(i <- 0 until ary.length)
        {
            val m=ssum(ary.length-i)+ssum(ary.length)-ssum(i)
            r(i) = 2*r(i)/m
        }
        return r
    }
    //mpm
    //maximam of the posterior marginals
    def mpm(ary:Array[Double], sampleRate:Double): Double =
    {
        val nary = nsdf(ary)

        var fl=false
        var first=true
        var localmax= -10.0
        var maxfreq= -1
        var allmax= -10.0
        var list=Vector[(Double, Int)]()
        for(i <- 0 until nary.length)
        {
            if(nary(i)<0)
            {
                first=false
                if(fl)
                {
                    fl=false
                    //push
                    list = list :+ ((localmax, maxfreq))
                    localmax= -10
                }
            }
            else
            {
                if(!first)
                {
                    fl=true
                    if(localmax<nary(i))
                    {
                        localmax=nary(i)
                        maxfreq=i
                    }
                    allmax=Math.max(allmax, nary(i))
                }
            }
        }
        //if(allmax<0.3) return 0.0
        for(i <- 0 until list.length)
        {
            val (f,s)=list(i)
            if(allmax*0.8<f)
            {
                return sampleRate/s
            }
        }
        return 0.0
    }
}
