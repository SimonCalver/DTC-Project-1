using UnityEngine;
using System.Collections;

public class WaterDetector : MonoBehaviour {

	public int index;
	public string type;

    void OnTriggerEnter2D(Collider2D Hit)
    {
        if (Hit.GetComponent<Rigidbody2D>() != null)
        {
			if(type == "WaterTrigger")
			{
				transform.parent.GetComponent<WorldManager> ().Splash (Mathf.Sqrt (Mathf.Pow (Hit.GetComponent<Rigidbody2D> ().velocity.x, 2) + Mathf.Pow (Hit.GetComponent<Rigidbody2D> ().velocity.y, 2)) * Hit.GetComponent<Rigidbody2D> ().mass / 40f, 
					index, Mathf.Atan2(Hit.GetComponent<Rigidbody2D> ().position.y,Hit.GetComponent<Rigidbody2D> ().position.x));
			}
			if(type == "LandTrigger")
			{
		//		transform.parent.GetComponent<LandCircle>().Splash(Hit.GetComponent<Rigidbody2D>().velocity.x*Hit.GetComponent<Rigidbody2D>().mass / 40f, Hit.GetComponent<Rigidbody2D>().velocity.y*Hit.GetComponent<Rigidbody2D>().mass / 40f,index);
		//		transform.parent.GetComponent<LandCircle>().impactCount = 5;
			}
		}
    }

    /*void OnTriggerStay2D(Collider2D Hit)
    {
        //print(Hit.name);
        if (Hit.rigidbody2D != null)
        {
            int points = Mathf.RoundToInt(Hit.transform.localScale.x * 15f);
            for (int i = 0; i < points; i++)
            {
                transform.parent.GetComponent<Water>().Splish(Hit.transform.position.x - Hit.transform.localScale.x + i * 2 * Hit.transform.localScale.x / points, Hit.rigidbody2D.mass * Hit.rigidbody2D.velocity.x / 10f / points * 2f);
            }
        }
    }*/

}
