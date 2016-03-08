using UnityEngine;
using System.Collections;

public class Clickspawn : MonoBehaviour {

    public GameObject Brick;
	
	void Update () {
        if (Input.GetMouseButtonDown(0))
        {
			Physics2D.gravity = new Vector3(0f,0.0f,0f);
            GameObject brick = Instantiate(Brick, Camera.main.ScreenToWorldPoint(Input.mousePosition + new Vector3(0,0,10)), Brick.transform.rotation) as GameObject;
            brick.transform.Rotate(0, 0, UnityEngine.Random.Range(0, 0));
			brick.transform.localScale = new Vector3(UnityEngine.Random.Range(0.5f, 7.5f),UnityEngine.Random.Range(0.5f, 7.5f),1)*0.4f;
            brick.GetComponent<Rigidbody2D>().mass = brick.transform.localScale.x * brick.transform.localScale.y*5f;
			float magsqr;
			Vector3 offset;
			offset = -brick.transform.position;
			offset.z = 0;
			magsqr = offset.sqrMagnitude;
			brick.GetComponent<Rigidbody2D>().AddForce((9.8f * offset) * 2.0f * brick.GetComponent<Rigidbody2D>().mass);
			Destroy(brick, 2.4f);

        }
	}
}
